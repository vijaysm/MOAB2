/*
 * Library for determining the type of a geometric model file
 *
 * Copyright 2006, Jason Kraftcheck (kraftche@cae.wisc.edu)
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#ifndef CGM2MOAB_HPP
#define CGM2MOAB_HPP

bool cgm2moab(MBInterface* iface,
              double dist_tol = 0.001,
              int norm_tol = 5,
              double len_tol = 0.0,
              int actuate_attribs = true);

#endif
