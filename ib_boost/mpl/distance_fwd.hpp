
#ifndef BOOST_MPL_DISTANCE_FWD_HPP_INCLUDED
#define BOOST_MPL_DISTANCE_FWD_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2000-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: distance_fwd.hpp 49267 2008-10-11 06:19:02Z agurtovoy $
// $Date: 2008-10-11 02:19:02 -0400 (Sat, 11 Oct 2008) $
// $Revision: 49267 $

#include <ib_boost/mpl/aux_/common_name_wknd.hpp>

namespace boost { namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(distance)

template< typename Tag > struct distance_impl;
template< typename First, typename Last > struct distance;

}}

#endif // BOOST_MPL_DISTANCE_FWD_HPP_INCLUDED
