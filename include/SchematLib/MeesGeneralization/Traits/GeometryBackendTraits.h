/*
  Copyright (C) 2018-2019 HERE Global B.V. and its affiliate(s).
All rights reserved.

This software and other materials contain proprietary information
controlled by HERE and are protected by applicable copyright legislation.
Any use and utilization of this software and other materials and
disclosure to any third parties is conditional upon having a separate
agreement with HERE for the access, use, utilization or disclosure of this
software. In the absence of such agreement, the use of the software is not
allowed.
 */
//
// Created by Mitra, Aniket on 29/01/2019.
//
#ifndef TULIB_TULIBGEOMETRY_BOOSTBACKEND_H
#define TULIB_TULIBGEOMETRY_BOOSTBACKEND_H

#include "tulib/logging.h"
#if CGAL_BACKEND_ENABLED
#include "tulib/geom/CGALTraits.h"
#else
#include "tulib/geom/BoostGeometryTraits.h"
#endif
#include "tulib/geom/GeometryInterface.h"
#include "tulib/metric/Norm.h"


struct GeometryKernel {

    #if CGAL_BACKEND_ENABLED
            //==============================
            // Define the Number Type
            // For the CGAL backend,
            // One can choose from the
            // following number types
            typedef long double NT;
            //typedef CGAL::Mpzf NT;
            //typedef CGAL::Gmpfr NT;
            //typedef CGAL::Gmpq NT;
            //==============================

            //==============================
            // Define the Dimensions
            // It is possible to choose
            // a higher dimension
            const static  size_t dimensions = 2;
            //const size_t dimensions = 7;
            //==============================

            //==============================
            // Define the Geometry Backend
            typedef tulib_support::CGALTraits<NT, dimensions> CGAL_GeometryBackend;
            //Using the Geometry Backend define the Tulib Geometry Kernel
            typedef tulib_core::TulibGeometryKernel<
                    typename CGAL_GeometryBackend::Wrapper_CGAL_Geometry> TulibGeometryKernel;
    #else
        //==============================
        // Define the Number Type
        // For the Boost Backend
        typedef long double NT;
        //==============================

        //==============================
        // Define the Dimensions
        // It is possible to choose
        // a higher dimension
        // Note: Boost Geometry Adapter supports geometry in two
        // dimensions only
        const static size_t dimensions = 2;
        //==============================

        //==============================
        // Define the Geometry Backend
        typedef tulib_support::BoostGeometryTraits<long double,dimensions> Boost_Geometry_Backend;
        //Using the Geometry Backend define the Tulib Geometry Kernel
        typedef tulib_core::TulibGeometryKernel<typename Boost_Geometry_Backend::Wrapper_Boost_Geometry> TulibGeometryKernel;
        //==============================
    #endif

    typedef tulib_support::FiniteNorm<TulibGeometryKernel,2> Norm;
};

#endif //TULIB_TULIBGEOMETRY_BOOSTBACKEND_H
