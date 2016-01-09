/*********************************************
Utility Library Function
Copyright(c) 2005, S. D. Rajan
All rights reserved

Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_CLOCK_H__
#define __RAJAN_CLOCK_H__

#include <time.h>
#include <string>

// -------------------------------------------
// Note all elapsed times are WALL CLOCK times
// -------------------------------------------

class CClock
{
    public:
        CClock ();				// ctor
        CClock (const CClock&); // overloaded ctor
        ~CClock();				// dtor
        
        // helper functions
        double DiffTime (const CClock&) const;
        double DiffTime () const;
        void GetDateTime (std::string& szDateTime) const;
        
        // modifier function
        void MarkTime ();

        // accessor function
        long GetMarkedTime () const;  

    private:
        time_t m_Time; // stores marked time since UTC 1/1/70
};

#endif
