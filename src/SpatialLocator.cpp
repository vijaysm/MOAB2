#include "moab/SpatialLocator.hpp"
#include "moab/Interface.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/BVHTree.hpp"

bool debug = true;

namespace moab 
{

    SpatialLocator::SpatialLocator(Interface *impl, Range &elems, Tree *tree, ElemEvaluator *eval) 
            : mbImpl(impl), myElems(elems), myDim(-1), myTree(tree), elemEval(eval), iCreatedTree(false)
    {
      create_tree();
      
      if (!elems.empty()) {
        myDim = mbImpl->dimension_from_handle(*elems.rbegin());
        ErrorCode rval = myTree->build_tree(myElems);
        if (MB_SUCCESS != rval) throw rval;
      }
    }

    void SpatialLocator::create_tree() 
    {
      if (myTree) return;
      
      if (myElems.empty() || mbImpl->type_from_handle(*myElems.rbegin()) == MBVERTEX) 
          // create a kdtree if only vertices
        myTree = new AdaptiveKDTree(mbImpl);
      else
          // otherwise a BVHtree, since it performs better for elements
        myTree = new BVHTree(mbImpl);

      iCreatedTree = true;
    }

    ErrorCode SpatialLocator::add_elems(Range &elems) 
    {
      if (elems.empty() ||
          mbImpl->dimension_from_handle(*elems.begin()) != mbImpl->dimension_from_handle(*elems.rbegin()))
        return MB_FAILURE;
  
      myDim = mbImpl->dimension_from_handle(*elems.begin());
      myElems = elems;

      ErrorCode rval = myTree->build_tree(myElems);
      return rval;
    }
    
#ifdef USE_MPI
    ErrorCode SpatialLocator::par_locate_points(Range &/*vertices*/,
                                                const double /*rel_iter_tol*/, const double /*abs_iter_tol*/,
                                                const double /*inside_tol*/)
    {
      return MB_UNSUPPORTED_OPERATION;
    }

    ErrorCode SpatialLocator::par_locate_points(const double */*pos*/, int /*num_points*/,
                                                const double /*rel_iter_tol*/, const double /*abs_iter_tol*/,
                                                const double /*inside_tol*/)
    {
      return MB_UNSUPPORTED_OPERATION;
    }
#endif
      
    ErrorCode SpatialLocator::locate_points(Range &verts,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol) 
    {
      assert(!verts.empty() && mbImpl->type_from_handle(*verts.rbegin()) == MBVERTEX);
      std::vector<double> pos(3*verts.size());
      ErrorCode rval = mbImpl->get_coords(verts, &pos[0]);
      if (MB_SUCCESS != rval) return rval;
      rval = locate_points(&pos[0], verts.size(), rel_iter_tol, abs_iter_tol, inside_tol);
      if (MB_SUCCESS != rval) return rval;
      
      return MB_SUCCESS;
    }
    
    ErrorCode SpatialLocator::locate_points(const double *pos, int num_points,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol) 
    {
        // initialize to tuple structure (p_ui, hs_ul, r[3]_d) (see header comments for locTable)
      locTable.initialize(1, 0, 1, 3, num_points);
      locTable.enableWriteAccess();

        // pass storage directly into locate_points, since we know those arrays are contiguous
      ErrorCode rval = locate_points(pos, num_points, (EntityHandle*)locTable.vul_wr, locTable.vr_wr, NULL, rel_iter_tol, abs_iter_tol,
                                     inside_tol);
      std::fill(locTable.vi_wr, locTable.vi_wr+num_points, 0);
      locTable.set_n(num_points);
      if (MB_SUCCESS != rval) return rval;
      
      return MB_SUCCESS;
    }
      
    ErrorCode SpatialLocator::locate_points(Range &verts,
                                            EntityHandle *ents, double *params, bool *is_inside,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol)
    {
      assert(!verts.empty() && mbImpl->type_from_handle(*verts.rbegin()) == MBVERTEX);
      std::vector<double> pos(3*verts.size());
      ErrorCode rval = mbImpl->get_coords(verts, &pos[0]);
      if (MB_SUCCESS != rval) return rval;
      return locate_points(&pos[0], verts.size(), ents, params, is_inside, rel_iter_tol, abs_iter_tol, inside_tol);
    }

    ErrorCode SpatialLocator::locate_points(const double *pos, int num_points,
                                            EntityHandle *ents, double *params, bool *is_inside,
                                            const double rel_iter_tol, const double abs_iter_tol, 
                                            const double inside_tol)
    {
      double tmp_abs_iter_tol = abs_iter_tol;
      if (rel_iter_tol && !tmp_abs_iter_tol) {
          // relative epsilon given, translate to absolute epsilon using box dimensions
        BoundBox box;
        myTree->get_bounding_box(box);
        tmp_abs_iter_tol = rel_iter_tol * box.diagonal_length();
      }
  
      EntityHandle closest_leaf;
      std::vector<double> dists;
      std::vector<EntityHandle> leaves;
      ErrorCode rval = MB_SUCCESS;

      for (int i = 0; i < num_points; i++) {
        int i3 = 3*i;
        ents[i] = 0;
        if (tmp_abs_iter_tol) {
          rval = myTree->distance_search(pos+i3, tmp_abs_iter_tol, leaves, tmp_abs_iter_tol, inside_tol, &dists);
          if (MB_SUCCESS != rval) return rval;
          if (!leaves.empty()) {
              // get closest leaf
            double min_dist = *dists.begin();
            closest_leaf = *leaves.begin();
            std::vector<EntityHandle>::iterator vit = leaves.begin()+1;
            std::vector<double>::iterator dit = dists.begin()+1;
            for (; vit != leaves.end() && min_dist; vit++, dit++) {
              if (*dit < min_dist) {
                min_dist = *dit;
                closest_leaf = *vit;
              }
            }
            dists.clear();
            leaves.clear();
          }
        }
        else {
          rval = myTree->point_search(pos+i3, closest_leaf);
          if (MB_ENTITY_NOT_FOUND == rval) closest_leaf = 0;
          else if (MB_SUCCESS != rval) return rval;
        }

          // if no ElemEvaluator, just return the box
        if (!elemEval) {
          ents[i] = closest_leaf;
          params[i3] = params[i3+1] = params[i3+2] = -2;
          if (is_inside && closest_leaf) is_inside[i] = true;
          continue;
        }
    
          // find natural coordinates of point in element(s) in that leaf
        CartVect tmp_nat_coords; 
        Range range_leaf;
        rval = mbImpl->get_entities_by_dimension(closest_leaf, myDim, range_leaf, false);
        if(rval != MB_SUCCESS) return rval;

          // loop over the range_leaf
        bool tmp_inside;
        bool *is_ptr = (is_inside ? is_inside+i : &tmp_inside);      
        *is_ptr = false;
        EntityHandle ent = 0;
        for(Range::iterator rit = range_leaf.begin(); rit != range_leaf.end(); rit++)
        {
          rval = elemEval->set_ent_handle(*rit); 
          if (MB_SUCCESS != rval) return rval;
          rval = elemEval->reverse_eval(pos+i3, tmp_abs_iter_tol, inside_tol, params+i3, is_ptr);
          if (MB_SUCCESS != rval) return rval;
          if (*is_ptr) {
            ent = *rit;
            break;
          }
        }
        if (debug && !ent) {
          std::cout << "Point " << i << " not found; point: (" 
                    << pos[i3] << "," << pos[i3+1] << "," << pos[i3+2] << ")" << std::endl;
          std::cout << "Source element candidates: " << std::endl;
          range_leaf.print("   ");
          for(Range::iterator rit = range_leaf.begin(); rit != range_leaf.end(); rit++)
          {
            std::cout << "Candidate " << CN::EntityTypeName(mbImpl->type_from_handle(*rit)) << " " << mbImpl->id_from_handle(*rit) << ": ";
            rval = elemEval->set_ent_handle(*rit); 
            if (MB_SUCCESS != rval) return rval;
            rval = elemEval->reverse_eval(pos+i3, tmp_abs_iter_tol, inside_tol, params+i3, is_ptr);
            if (MB_SUCCESS != rval) return rval;
            std::cout << "Parameters: (" << params[i3] << "," << params[i3+1] << "," << params[i3+2] << ")" 
                      << " inside = " << *is_ptr << std::endl;
          }
        }
        ents[i] = ent;
      }

      return MB_SUCCESS;
    }
    
        /* Count the number of located points in locTable
         * Return the number of entries in locTable that have non-zero entity handles, which
         * represents the number of points in targetEnts that were inside one element in sourceEnts
         *
         */
    int SpatialLocator::local_num_located() 
    {
      int num_located = locTable.get_n() - std::count(locTable.vul_rd, locTable.vul_rd+locTable.get_n(), 0);
      if (num_located != (int)locTable.get_n()) {
        unsigned long *nl = std::find(locTable.vul_rd, locTable.vul_rd+locTable.get_n(), 0);
        if (nl) {
          int idx = nl - locTable.vul_rd;
          if (idx) {}
        }
      }
      return num_located;
    }

        /* Count the number of located points in parLocTable
         * Return the number of entries in parLocTable that have a non-negative index in on a remote
         * proc in parLocTable, which gives the number of points located in at least one element in a
         * remote proc's sourceEnts.
         */
    int SpatialLocator::remote_num_located()
    {
      int located = 0;
      for (unsigned int i = 0; i < parLocTable.get_n(); i++)
        if (parLocTable.vi_rd[2*i] != -1) located++;
      return located;
    }
} // namespace moab

