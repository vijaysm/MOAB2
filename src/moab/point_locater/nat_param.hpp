#ifndef NAT_PARAM_HPP
#define NAT_PARAM_HPP
template< typename Entities, typename Coordinates>
ErrorCode nat_param( double query_point[3], 
                     Entities & entities, 
                     Coordinates & nat_coords,
                     double epsilon = 0.0){
  AdaptiveKDTreeIter treeiter;
  ErrorCode result = tree->get_tree_iterator(local_root, treeiter); 
  if (MB_SUCCESS != result) {
    std::cout << "Problems getting iterator" << std::endl;
    return result;
  }

  EntityHandle closest_leaf;
  if (epsilon) {
    Vector_double dists;
    Entities leaves;
    result = tree->leaves_within_distance( local_root, query_points, 
					   epsilon, leaves, &dists);
    if (leaves.empty()){ 
      // not found returns success here, with empty list, just like case with no epsilon
      return MB_SUCCESS;
      }
      // get closest leaf
    //TODO: replace all this with a call to std::min() and std::dist()
    double min_dist = *dists.begin();
    closest_leaf = *leaves.begin();
    Entities::iterator vit = leaves.begin()+1;
    Vector_double::iterator dit = dists.begin()+1;
    for (; vit != leaves.end() && min_dist; vit++, dit++) {
      if (*dit < min_dist) {
        min_dist = *dit;
        closest_leaf = *vit;
      }
    }
  }
  else {
    switch( result=tree->leaf_containing_point( local_root, 
						query_points, 
						treeiter)){
    //point is outside of tree's bounding box
    case: MB_ENTITY_NOT_FOUND 
      return MB_SUCCESS;
    case: MB_SUCCESS 
      closest_leaf = treeiter.handle();
      break;
    default:
	return result;
   }
  }

  // find natural coordinates of point in element(s) in that leaf
  CartVect tmp_nat_coords; 
  Range range_leaf;
  if (!impl.get_entities_by_dimension(closest_leaf, 3, range_leaf, false)){
  	std::cout << "Problem getting leaf in a range" << std::endl;
   }

  // loop over the range_leaf
  for( Range::iterator iter = range_leaf.begin(); 
		       iter != range_leaf.end(); 
					 iter++)
  {
    //test to find out in which entity the point is
    //get the EntityType and create the appropriate Element::Map subtype
    // if spectral, do not need coordinates, just the GL points
    EntityType etype = impl.type_from_handle(*iter);
    if (NULL!= this->_spectralSource && etype==MBHEX)
    {
      EntityHandle eh = *iter;
      const double * xval;
      const double * yval;
      const double * zval;
      ErrorCode rval = impl.tag_get_by_ptr(_xm1Tag, &eh, 1,(const void **) &xval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get xm1 values \n";
        return MB_FAILURE;
      }
      rval = impl.tag_get_by_ptr(_ym1Tag, &eh, 1, (const void **)&yval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get ym1 values \n";
        return MB_FAILURE;
      }
      rval = impl.tag_get_by_ptr(_zm1Tag, &eh, 1, (const void **)&zval );
      if (moab::MB_SUCCESS != rval)
      {
        std::cout << "can't get zm1 values \n";
        return MB_FAILURE;
      }
      Element::SpectralHex * spcHex = ( Element::SpectralHex * ) _spectralSource;

      spcHex->set_gl_points((double*)xval, (double*)yval, (double*)zval);
      try{
        tmp_nat_coords =spcHex->ievaluate(CartVect(query_points));
      }
      catch (Element::Map::EvaluationError) {
        std::cout << "point "<< query_points[0] << " " << query_points[1] << " " << query_points[2] <<
            " is not converging inside hex " << impl.id_from_handle(eh) << "\n";
        continue; // it is possible that the point is outside, so it will not converge
      }

    }
    else
    {
      const EntityHandle *connect;
      int num_connect;

        //get connectivity
      result = impl.get_connectivity(*iter, connect, num_connect, true);

        //get coordinates of the vertices
      std::vector<CartVect> coords_vert( num_connect);
      result = impl.get_coords( connect, num_connect, &(coords_vert[0][0]));
      if (MB_SUCCESS != result) {
        std::cout << "Problems getting coordinates of vertices\n";
        return result;
      }

      if (etype == MBHEX) {
        Element::LinearHex hexmap(coords_vert);
        try {
          tmp_nat_coords = hexmap.ievaluate(CartVect(query_points), epsilon);
        }
        catch (Element::Map::EvaluationError) {
          continue;
        }
      }
      else if (etype == MBTET){
        Element::LinearTet tetmap(coords_vert);
        try {
          tmp_nat_coords = tetmap.ievaluate(CartVect(query_points));
        }
        catch (Element::Map::EvaluationError) {
          continue;
        }
      }
      else {
        std::cout << "Entity not Hex or Tet" << std::endl;
        continue;
      }
    }
      //if we get here then we've found the coordinates.
      //save them and the entity and return success.
    entities.push_back(*iter);
    nat_coords.push_back(tmp_nat_coords);
    return MB_SUCCESS;
  }

  //didn't find any elements containing the point
  return MB_SUCCESS;
}
#endif //NAT_PARAM_HPP
