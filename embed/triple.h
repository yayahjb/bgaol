/*
 *  triple.h
 *  
 *
 *  Created by Herbert J. Bernstein on 10/26/10.
 *  Copyright 2010 Herbert J. Bernstein
 *
 */
 
/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CVector API UNDER THE TERMS OF THE LGPL   *
 *                                                                    *
 **********************************************************************/

/************************* LGPL NOTICES *******************************
 *                                                                    *
 * This library is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU Lesser General Public         *
 * License as published by the Free Software Foundation; either       *
 * version 2.1 of the License, or (at your option) any later version. *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
 * Lesser General Public License for more details.                    *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License along with this library; if not, write to the Free         *
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
 * MA  02110-1301  USA                                                *
 *                                                                    *
 **********************************************************************/

/************************* PAIR.H NOTICE ******************************
 * triple.h was derived from pair.h                                   *
 *  which was subject to the following notice:                        *
 *                                                                    *
 * Copyright (c) 1994                                                 *
 * Hewlett-Packard Company                                            *
 *                                                                    *
 * Permission to use, copy, modify, distribute and sell this software *
 * and its documentation for any purpose is hereby granted without    *
 * fee, provided that the above copyright notice appear in all copies *
 * and that both that copyright notice and this permission notice     *
 * appear in supporting documentation.  Hewlett-Packard Company makes *
 * no representations about the suitability of this software for any  *
 * purpose.  It is provided "as is" without express or implied        *
 * warranty.                                                          *
 *                                                                    *
 **********************************************************************/

#ifndef TRIPLE_H
#define TRIPLE_H

template <class TR1, class TR2, class TR3>
class triple {
private:
    TR1 m_first;
    TR2 m_second;
    TR3 m_third;
public:
    triple() : m_first(TR1()), m_second(TR2()), m_third(TR3()) {}
    triple(const TR1& first, const TR2& second, const TR3& third ) : m_first(first), m_second(second), m_third(third) {}
    inline TR1 GetFirst( void ) const  { return m_first;  }  
    inline TR2 GetSecond( void ) const { return m_second; }  
    inline TR3 GetThird( void ) const  { return m_third;  }  

inline bool operator==(const triple<TR1, TR2, TR3>& rhs) { 
    return m_first == rhs.GetFirst() && m_second == rhs.m_second && m_third == rhs.m_third; 
}

inline bool operator<(const triple<TR1, TR2, TR3>& rhs) { 
    return m_first < rhs.GetFirst() 
      || (m_first == rhs.GetFirst() && m_second < rhs.GetSecond)
      || (m_first == rhs.GetFirst() && m_second < rhs.GetSecond && m_third < rhs.GetThird); 
}

};

template <class TR1, class TR2, class TR3>
inline triple<TR1, TR2, TR3> make_triple(const TR1& first, const TR2& second, const TR3& third) {
    return triple<TR1, TR2, TR3>(first, second, third);
}


#endif
