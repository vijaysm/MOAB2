#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/CartVect.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/CN.hpp"
#include "ElemUtil.hpp"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>
#include <string>

// Different platforms follow different conventions for usage
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

using namespace moab;
#define CHK(r,s) do { if (MB_SUCCESS != (r)) {fail( (r), (s), __FILE__, __LINE__ ); return (r);}} while(false)

const int ELEMEVAL = 0, ELEMUTIL = 1;

static void fail( ErrorCode error_code, const char *str, const char* file_name, int line_number )
{
  std::cerr << str << ", line " << line_number << " of file " << file_name << ", error code " << error_code << std::endl;
}

double mytime();

ErrorCode get_ents(Interface &mbi, std::string &filename, int &dim, Range &elems, EntityType &tp, int &nv) 
{
  ErrorCode rval = mbi.load_file(filename.c_str(), 0);
  while (elems.empty() && dim >= 1) {
    rval = mbi.get_entities_by_dimension(0, dim, elems); CHK(rval, "get_entities_by_dimension");
    if (elems.empty()) dim--;
  }
  if (elems.empty()) {
    CHK(MB_FAILURE, "No elements in file");
  }

    // check to see they're all the same type & #vertices
  tp = mbi.type_from_handle(*elems.begin());
  EntityType tp2 = mbi.type_from_handle(*elems.rbegin());
  if (tp != tp2) 
    CHK(MB_FAILURE, "Elements must have same type");

  int nv2;
  const EntityHandle *c1;
  rval = mbi.get_connectivity(*elems.begin(), c1, nv); CHK(rval, "get_connectivity");
  rval = mbi.get_connectivity(*elems.rbegin(), c1, nv2); CHK(rval, "get_connectivity");
  if (nv != nv2) {
    CHK(MB_FAILURE, "Elements must have same #vertices");
  }
        
  return MB_SUCCESS;
}

void parse_options(ProgOptions &opts, int &dim, std::string &filename)
{
  opts.addOpt<int>(std::string("dim"), std::string("Evaluate 1d/2d/3d elements (default: maximal dimension in mesh)"),
                   &dim, ProgOptions::int_flag);
  opts.addOpt<std::string>(std::string("filename,f"), std::string("Filename containing mesh"), &filename);
}

ErrorCode get_elem_map(EntityType tp, std::vector<CartVect> &vcoords, int nconn, Element::Map *&elemmap) 
{
  switch (tp) {
    case MBHEX:
        if (nconn == 8) {
          elemmap = new Element::LinearHex(vcoords);
          break;
        }
        else if (nconn == 27) {
          elemmap = new Element::QuadraticHex(vcoords);
          break;
        }
        else return MB_FAILURE;
    case MBTET:
        if (nconn == 4) {
          elemmap = new Element::LinearTet(vcoords);
          break;
        }
        else return MB_FAILURE;
    case MBQUAD:
        if (nconn == 4) {
          elemmap = new Element::LinearQuad(vcoords);
          break;
        }
        else return MB_FAILURE;
    default:
        return MB_FAILURE;
  }
  
  return MB_SUCCESS;
}

ErrorCode time_forward_eval(Interface *mbi, int method, Range &elems, 
                            std::vector<CartVect> &params, std::vector<CartVect> &coords, 
                            double &evtime) 
{
  evtime = mytime();
  ErrorCode rval;
  Range::iterator rit;
  unsigned int i;
  if (ELEMEVAL == method) {
      // construct ElemEvaluator
    EvalSet eset;
    ElemEvaluator eeval(mbi);
    eeval.set_eval_set(*elems.begin());
    eeval.set_tag_handle(0, 0); // indicates coordinates as the field to evaluate
  
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      eeval.set_ent_handle(*rit);
      rval = eeval.eval(params[i].array(), coords[i].array(), 3);
#ifndef NDEBUG
      if (MB_SUCCESS != rval) return rval;
#endif
    }
  }
  else if (ELEMUTIL == method) {
    std::vector<CartVect> vcoords(CN::MAX_NODES_PER_ELEMENT);
    const EntityHandle *connect;
    int nconn;
    Element::Map *elemmap = NULL;
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      rval = mbi->get_connectivity(*rit, connect, nconn); CHK(rval, "get_connectivity");
      rval = mbi->get_coords(connect, nconn, vcoords[0].array()); CHK(rval, "get_coords");
      rval = get_elem_map(mbi->type_from_handle(*rit), vcoords, nconn, elemmap); CHK(rval, "get_elem_map");
      coords[i] = elemmap->evaluate(params[i]);
    }
  }
  
  evtime = mytime() - evtime;
  return MB_SUCCESS;
}

ErrorCode time_reverse_eval(Interface *mbi, int method, Range &elems, 
                            std::vector<CartVect> &coords, std::vector<CartVect> &params, 
                            double &retime) 
{
  retime = mytime();
  ErrorCode rval;
  Range::iterator rit;
  unsigned int i;
  if (ELEMEVAL == method) {
    EvalSet eset;
    ElemEvaluator eeval(mbi);
    eeval.set_eval_set(*elems.begin());
    eeval.set_tag_handle(0, 0); // indicates coordinates as the field to evaluate    
    int ins;
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      eeval.set_ent_handle(*rit);
      rval = eeval.reverse_eval(coords[i].array(), 1.0e-10, 1.0e-6, params[i].array(), &ins);
      assert(ins);
#ifndef NDEBUG
      if (MB_SUCCESS != rval) return rval;
#endif
    }
  }
  else if (ELEMUTIL == method) {
    std::vector<CartVect> vcoords(CN::MAX_NODES_PER_ELEMENT);
    const EntityHandle *connect;
    int nconn;
    Element::Map *elemmap = NULL;
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      rval = mbi->get_connectivity(*rit, connect, nconn); CHK(rval, "get_connectivity");
      rval = mbi->get_coords(connect, nconn, vcoords[0].array()); CHK(rval, "get_coords");
      rval = get_elem_map(mbi->type_from_handle(*rit), vcoords, nconn, elemmap); CHK(rval, "get_elem_map");
      coords[i] = elemmap->ievaluate(coords[i], 1.0e-6);
    }
  }
  retime = mytime() - retime;
  return MB_SUCCESS;
}

ErrorCode time_jacobian(Interface *mbi, int method, Range &elems, std::vector<CartVect> &params, 
                        double &jactime) 
{
  jactime = mytime();
  ErrorCode rval;
  Range::iterator rit;
  unsigned int i;
  Matrix3 jac;
  if (ELEMEVAL == method) {
    EvalSet eset;
    ElemEvaluator eeval(mbi);
    eeval.set_eval_set(*elems.begin());
    eeval.set_tag_handle(0, 0); // indicates coordinates as the field to evaluate    
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      eeval.set_ent_handle(*rit);
      rval = eeval.jacobian(params[i].array(), jac.array());
#ifndef NDEBUG
      if (MB_SUCCESS != rval) return rval;
#endif
    }
  }
  else if (ELEMUTIL == method) {
    std::vector<CartVect> vcoords(CN::MAX_NODES_PER_ELEMENT);
    const EntityHandle *connect;
    int nconn;
    Element::Map *elemmap = NULL;
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      rval = mbi->get_connectivity(*rit, connect, nconn); CHK(rval, "get_connectivity");
      rval = mbi->get_coords(connect, nconn, vcoords[0].array()); CHK(rval, "get_coords");
      rval = get_elem_map(mbi->type_from_handle(*rit), vcoords, nconn, elemmap); CHK(rval, "get_elem_map");
      jac = elemmap->jacobian(params[i]);
    }
  }
  jactime = mytime() - jactime;
  return MB_SUCCESS;
}

ErrorCode time_integrate(Interface *mbi, int method, Tag tag, Range &elems, double &inttime) 
{
  inttime = mytime();
  ErrorCode rval;
  double integral;
  if (ELEMEVAL == method) {
    EvalSet eset;
    ElemEvaluator eeval(mbi);
    eeval.set_eval_set(*elems.begin());
    eeval.set_tag_handle(0, 0); // indicates coordinates as the field to evaluate    
    rval = eeval.set_tag_handle(tag, 0); CHK(rval, "set_tag_handle");
    for (Range::iterator rit = elems.begin(); rit != elems.end(); ++rit) {
      eeval.set_ent_handle(*rit);
      rval = eeval.integrate(&integral);
#ifndef NDEBUG
      if (MB_SUCCESS != rval) return rval;
#endif
    }
  }
  else if (ELEMUTIL == method) {
    std::vector<CartVect> vcoords(CN::MAX_NODES_PER_ELEMENT);
    std::vector<double> tagval(CN::MAX_NODES_PER_ELEMENT);
    const EntityHandle *connect;
    int nconn;
    Element::Map *elemmap = NULL;
    Range::iterator rit;
    unsigned int i;
    for (rit = elems.begin(), i = 0; rit != elems.end(); ++rit, i++) {
      rval = mbi->get_connectivity(*rit, connect, nconn); CHK(rval, "get_connectivity");
      rval = mbi->get_coords(connect, nconn, vcoords[0].array()); CHK(rval, "get_coords");
      rval = get_elem_map(mbi->type_from_handle(*rit), vcoords, nconn, elemmap); CHK(rval, "get_elem_map");
      rval = mbi->tag_get_data(tag, connect, nconn, &tagval[0]); CHK(rval, "tag_get_data");
      integral = elemmap->integrate_scalar_field(&tagval[0]);
    }
  }
  inttime = mytime() - inttime;
  std::cout << "integral = " << integral << std::endl;
  return MB_SUCCESS;
}

ErrorCode put_random_field(Interface &mbi, Tag &tag, Range &elems) 
{
  Range verts;
  ErrorCode rval = mbi.get_adjacencies(elems, 0, false, verts, Interface::UNION); CHK(rval, "get_adjacencies");
  rval = mbi.tag_get_handle("mytag", 1, MB_TYPE_DOUBLE, tag, MB_TAG_CREAT | MB_TAG_DENSE); CHK(rval, "tag_get_handle");
  std::vector<double> tag_vals(verts.size());
  for (unsigned int i = 0; i < verts.size(); i++)
    tag_vals[i] = ((double)rand())/RAND_MAX;
  rval = mbi.tag_set_data(tag, verts, &tag_vals[0]);
  return rval;
}

ErrorCode elem_evals(Interface *mbi, int method, Range &elems, Tag tag,
                     std::vector<CartVect> &params, std::vector<CartVect> &coords, 
                     double &evtime, double &retime, double &jactime, double &inttime) 
{
  evtime = 0, retime = 0, jactime = 0, inttime = 0; // initializations to avoid compiler warning
  
    // time/test forward evaluation, putting results into vector
  ErrorCode rval = time_forward_eval(mbi, method, elems, params, coords, evtime); CHK(rval, "time_forward_eval");

    // time/test reverse evaluation, putting results into vector
  rval = time_reverse_eval(mbi, method, elems, coords, params, retime); CHK(rval, "time_reverse_eval");

    // time/test Jacobian evaluation
  rval = time_jacobian(mbi, method, elems, params, jactime); CHK(rval, "time_jacobian");

    // time/test integration
  rval = time_integrate(mbi, method, tag, elems, inttime); CHK(rval, "time_integrate");

  return rval;
}

int main( int argc, char* argv[] )
{
    // parse options
  ProgOptions opts;
  std::string filename;
  int dim = 3;
  parse_options(opts, dim, filename);
  opts.parseCommandLine(argc, argv);
  if (filename.empty()) CHK(MB_FAILURE, "No file specified");
  else if (dim < 1 || dim > 3) CHK(MB_FAILURE, "Dimension must be > 0 and <= 3");
  
    // read mesh file & gather element handles
  Core mbi;
  Range elems;
  int nv;
  EntityType tp;
  ErrorCode rval = get_ents(mbi, filename, dim, elems, tp, nv);
  if (MB_SUCCESS != rval) return 1;

    // construct (random) parameter values for queries
  unsigned int ne = elems.size();
  std::vector<CartVect> params(ne), coords(ne);
  srand(time(NULL));
  for (unsigned int i = 0; i < ne; i++) {
    params[i][0] = -1 + 2*((double)rand())/RAND_MAX;
    if (dim > 1) params[i][1] = -1 + 2*((double)rand())/RAND_MAX;
    if (dim > 2) params[i][2] = -1 + 2*((double)rand())/RAND_MAX;
  }

    // put random field on vertices
  Tag tag;
  rval = put_random_field(mbi, tag, elems);
  CHK(rval, "put_random_field");
  double evtime[2], retime[2], jactime[2], inttime[2]; // initializations to avoid compiler warning

  rval = elem_evals(&mbi, ELEMEVAL, elems, tag, params, coords,
                    evtime[0], retime[0], jactime[0], inttime[0]);
  CHK(rval, "new elem_evals");
  
  rval = elem_evals(&mbi, ELEMUTIL, elems, tag, params, coords, 
                    evtime[1], retime[1], jactime[1], inttime[1]);
  CHK(rval, "old elem_evals");
  
  std::cout << filename << ": " << elems.size() << " " << CN::EntityTypeName(tp)
            << " elements, " << nv << " vertices per element." << std::endl << std::endl;
  std::cout << "New, old element evaluation code:" << std::endl;
  std::cout << "Evaluation type, time, time per element:" << std::endl;
  std::cout << "                             New                   Old            (New/Old)*100" << std::endl;
  std::cout << "Forward evaluation " << evtime[0] << ", " << evtime[0] / elems.size()
            << "    " << evtime[1] << ", " << evtime[1] / elems.size() 
            << "    " << (evtime[0]/(evtime[1]?evtime[1]:1))*100.0 << std::endl;
  std::cout << "Reverse evaluation " << retime[0] << ", " << retime[0] / elems.size()
            << "    " << retime[1] << ", " << retime[1] / elems.size() 
            << "    " << (retime[0]/(retime[1]?retime[1]:1))*100.0 << std::endl;
  std::cout << "Jacobian           " << jactime[0] << ", " << jactime[0] / elems.size()
            << "    " << jactime[1] << ", " << jactime[1] / elems.size() 
            << "    " << (jactime[0]/(jactime[1]?jactime[1]:1))*100.0 << std::endl;
  std::cout << "Integration        " << inttime[0] << ", " << inttime[0] / elems.size()
            << "    " << inttime[1] << ", " << inttime[1] / elems.size() 
            << "    " << (inttime[0]/(inttime[1]?inttime[1]:1))*100.0 << std::endl;
}

#if defined(_MSC_VER) || defined(__MINGW32__)
double mytime2(double &tot_time, double &utime, double &stime, long &imem, long &rmem) 
{
  utime = (double)clock() / CLOCKS_PER_SEC;
  tot_time = stime = 0;
  imem = rmem = 0;
  return tot_time;
}
#else
double mytime2(double &tot_time, double &utime, double &stime, long &imem, long &rmem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  utime = (double)r_usage.ru_utime.tv_sec +
     ((double)r_usage.ru_utime.tv_usec/1.e6);
  stime = (double)r_usage.ru_stime.tv_sec +
     ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = utime + stime;
#ifndef LINUX
  imem = r_usage.ru_idrss;
  rmem = r_usage.ru_maxrss;
#else
  system("ps o args,drs,rss | grep perf | grep -v grep");  // RedHat 9.0 doesnt fill in actual memory data 
  imem = rmem = 0;
#endif
  return tot_time;
}
#endif
double mytime() 
{
  double ttime, utime, stime;
  long imem, rmem;
  return mytime2(ttime, utime, stime, imem, rmem);
}
