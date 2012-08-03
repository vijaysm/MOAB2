/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

//Point Locater
#include "moab/point_locater/point_locater.hpp"


//STL
#include <string>
#include <vector>
#include <sstream>

// default types.. whatevs.
typedef std::vector< int> Elements;
typedef std::vector< int> Tree;
typedef int Communicator; 
typedef moab::Point_search< Elements, Tree, Communicator> Point_search;

int main(int argc, char* argv[]){
	Point_search locator();
}
