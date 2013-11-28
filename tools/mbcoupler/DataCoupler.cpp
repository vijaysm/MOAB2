#include "DataCoupler.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/Tree.hpp"
#include "moab/TupleList.hpp"
#include "moab/SpatialLocator.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/Error.hpp"

#include "iostream"
#include <stdio.h>
#include <algorithm>
#include <sstream>

#include "assert.h"

namespace moab {

bool debug = false;

DataCoupler::DataCoupler(Interface *impl,
                         ParallelComm *pc,
                         Range &source_ents,
                         int coupler_id,
                         bool init_locator,
                         int dim)
        : mbImpl(impl), myPcomm(pc), myId(coupler_id), myDim(dim)
{
  assert(NULL != mbImpl && (myPcomm || !source_ents.empty()));

    // now initialize the tree
  if (init_locator) {
    myLocator = new SpatialLocator(mbImpl, source_ents);
    myLocator->elem_eval(new ElemEvaluator(mbImpl));

      // initialize element evaluator with the default for the entity types in source_ents; 
      // can be replaced later by application if desired
    if (!source_ents.empty()) {
      Range::pair_iterator pit = source_ents.pair_begin();
      EntityType last_type = MBMAXTYPE;
      for (; pit != source_ents.pair_end(); pit++) {
        EntityType this_type = mbImpl->type_from_handle(pit->first);
        if (last_type == this_type) continue;
        myLocator->elem_eval()->set_eval_set(pit->first);
        last_type = this_type;
      }
    }
  }
  
  if (-1 == dim && !source_ents.empty()) 
    dim = mbImpl->dimension_from_handle(*source_ents.rbegin());

  ErrorCode rval = impl->query_interface(mError);
  if (MB_SUCCESS != rval) throw(rval);
}

  /* destructor
   */
DataCoupler::~DataCoupler()
{
  delete myLocator;
}

ErrorCode DataCoupler::locate_points(Range &targ_ents,
                                     double rel_eps, 
                                     double abs_eps)
{
  targetEnts = targ_ents;
  
  return myLocator->locate_points(targ_ents, rel_eps, abs_eps);
}

ErrorCode DataCoupler::locate_points(double *xyz, int num_points,
                                 double rel_eps, 
                                 double abs_eps)
{
  return myLocator->locate_points(xyz, num_points, rel_eps, abs_eps);
}

ErrorCode DataCoupler::interpolate(/*DataCoupler::Method*/ int method,
                                   const std::string &interp_tag,
                                   double *interp_vals,
                                   std::vector<int> *point_indices,
                                   bool normalize)
{
    // tag name input, translate to tag handle and pass down the chain

    // not inlined because of call to set_last_error, class Error isn't in public interface
  Tag tag;
  ErrorCode result = mbImpl->tag_get_handle(interp_tag.c_str(), tag);
  if (MB_SUCCESS != result) {
    std::ostringstream str;
    str << "Failed to get handle for interpolation tag \"" << interp_tag << "\"";
    mError->set_last_error(str.str());
    return result;
  }
  return interpolate(method, tag, interp_vals, point_indices, normalize);
}
  
ErrorCode DataCoupler::interpolate(/*DataCoupler::Method*/ int *methods,
                                   Tag *tags,
                                   int *points_per_method,
                                   int num_methods,
                                   double *interp_vals,
                                   std::vector<int> *point_indices,
                                   bool /*normalize*/)
{
    // lowest-level interpolate function, does actual interpolation using calls to ElemEvaluator

  ErrorCode result = MB_SUCCESS;

  unsigned int pts_total = 0;
  for (int i = 0; i < num_methods; i++) pts_total += (points_per_method ? points_per_method[i] : targetEnts.size());

  unsigned int num_indices = (point_indices ? point_indices->size() : targetEnts.size());

  int max_tsize = -1;
  for (int i = 0; i < num_methods; i++) {
    int tmp_tsize;
    result = mbImpl->tag_get_length(tags[i], tmp_tsize);
    if (MB_SUCCESS != result) return MB_FAILURE;
    max_tsize = std::max(max_tsize, tmp_tsize);
  }

    // if tl was passed in non-NULL, just have those points, otherwise have targetPts plus
    // locally mapped pts
  if (pts_total != num_indices)
    return MB_FAILURE;

  if (myPcomm) {
      // TL to send interpolation indices to target procs
      // Tuple structure: (pto_i, ridx_i, lidx_i, meth_i, tagidx_i, interp_val[max_tsize]_d)
    TupleList tinterp;
    tinterp.initialize(5, 0, 0, max_tsize, num_indices);
    int t = 0;
    tinterp.enableWriteAccess();
    for (int i = 0; i < num_methods; i++) {
      int num_points = (points_per_method ? points_per_method[i] : targetEnts.size());
      for (int j = 0; j < num_points; j++) {
        int idx = (point_indices ? (*point_indices)[j] : j);

          // remote proc/idx from myLocator->parLocTable
        tinterp.vi_wr[5*t]   = myLocator->par_loc_table().vi_rd[2*idx]; // proc
        tinterp.vi_wr[5*t+1] = myLocator->par_loc_table().vi_rd[2*idx+1]; // remote idx
  
          // local entity index, tag/method index from my data
        tinterp.vi_wr[5*t+2] = idx;
        tinterp.vi_wr[5*t+3] = methods[i];
        tinterp.vi_wr[5*t+4] = i;
        tinterp.inc_n();
        t++;
      }
    }

      // scatter/gather interpolation points
    myPcomm->proc_config().crystal_router()->gs_transfer(1, tinterp, 0);

      // perform interpolation on local source mesh; put results into
      // tinterp.vr_wr

    for (unsigned int i = 0; i < tinterp.get_n(); i++) {
      int lidx = tinterp.vi_rd[5*i+1];
//    /*Method*/ int method = (/*Method*/ int)tinterp.vi_rd[5*i+3];
      Tag tag = tags[tinterp.vi_rd[5*i+4]];

      myLocator->elem_eval()->set_tag_handle(tag);
      myLocator->elem_eval()->set_ent_handle(myLocator->loc_table().vul_rd[lidx]);
      result = myLocator->elem_eval()->eval(myLocator->loc_table().vr_rd+3*lidx, tinterp.vr_rd+i*max_tsize);
      if (MB_SUCCESS != result) return result;
    }

      // scatter/gather interpolation data
    myPcomm->proc_config().crystal_router()->gs_transfer(1, tinterp, 0);

      // copy the interpolated field as a unit
    std::copy(tinterp.vr_rd, tinterp.vr_rd+tinterp.get_n()*max_tsize, interp_vals);
  }
  else {
    std::vector<double> tmp_vals;
    std::vector<EntityHandle> tmp_ents;
    double *tmp_dbl = interp_vals;
    for (int i = 0; i < num_methods; i++) {
      int num_points = (points_per_method ? points_per_method[i] : targetEnts.size());

        // interpolated data is tsize long, which is either max size (if data passed back to caller in tinterp)
        // or tag size (if data will be set on entities, in which case it shouldn't have padding)
      int tsize = max_tsize, tsize_bytes = 0;
      if (!interp_vals) {
        tmp_vals.resize(num_points*max_tsize);
        tmp_dbl = &tmp_vals[0];
        tmp_ents.resize(num_points);
        result = mbImpl->tag_get_length(tags[i], tsize);
        result = mbImpl->tag_get_bytes(tags[i], tsize_bytes);
      }
      
      for (int j = 0; j < num_points; j++) {
        int lidx;
        if (point_indices) {
          lidx = (*point_indices)[j];
        }
        else {
          lidx = j;
        }

        myLocator->elem_eval()->set_tag_handle(tags[i]);
        myLocator->elem_eval()->set_ent_handle(myLocator->loc_table().vul_rd[lidx]);
        if (!interp_vals) tmp_ents[j] = targetEnts[lidx]; // could be performance-sensitive, thus the if test
        result = myLocator->elem_eval()->eval(myLocator->loc_table().vr_rd+3*lidx, tmp_dbl);
        tmp_dbl += tsize;
        if (MB_SUCCESS != result) return result;
      } // for j

      if (!interp_vals) {
          // set tags on tmp_ents; data is already w/o padding, due to tsize setting above
        result = mbImpl->tag_set_data(tags[i], &tmp_ents[0], tmp_ents.size(), &tmp_vals[0]);
        if (MB_SUCCESS != result) return result;
      }

    } // for i
  } // if myPcomm
  
      // done
  return MB_SUCCESS;
}

} // namespace_moab
