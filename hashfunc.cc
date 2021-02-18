//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		CHashFunc: Class for Universal hash functions
//
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details (www.gnu.org).
//
//*****************************************************************/

#include "hashfunc.hh"
#include <cassert>
//#include "RandomLib/Random.hpp"

// constructor	
CHashFunc::CHashFunc(
					 unsigned int t, 
					 unsigned int n, 
					 unsigned int c) 
: _m1(0), _m2(0), _t(t), _n(n), _a1(NULL), _a2(NULL), _c(c)
{ 
	UHashfunc_init(t, n, c);
}

// destructor
CHashFunc::~CHashFunc()
{
	delete[] _a1;
	delete[] _a2;
}

void
CHashFunc::UHashfunc_init(
						  unsigned int t, 
						  unsigned int n, 
						  unsigned int c)
{	
	// Init member variables
	_t = t;
	_n = n;
	_c = c;
	
	// Get the prime number which is larger than t*n 
	unsigned long long top = _t*_n;
	
	unsigned long long p = 0;
	unsigned int mul = 1;
	do {		
		unsigned from = 100 * mul;
		p = GetPrime(top, from);
		++mul;
	} while (p == 0);
	
	_m1 = p;   	// t*n ~~ m1
	
	unsigned long long top2 = _c*_t*_n;
	unsigned long long p2 = 0;
	mul = 1;
	do {		
		unsigned from = 100 * mul;
		p2 = GetPrime(top2, from);
		++mul;
	} while (p2 == 0);
	
	_m2 = p2; // m2 > c*t*n --> to avoid double collision ==> I just use _c*top for _m2
    if(_a1 != NULL)
        delete[] _a1;
    if(_a2 != NULL)
        delete[] _a2;
	_a1 = new unsigned long long[_n];
	_a2 = new unsigned long long[_n];
	
	
	// generate n random numbers between 0 and m1-1
	// for hash function1 and hash function2
	// rand() % 48       
	// random number between 0~47
	//RandomLib::Random rnd;		// r created with random seed

    srand(1);
	for (unsigned int i=0; i<_n; ++i) {
		//_a1[i] = rnd.Integer<unsigned long long>(_m1-1);
		_a1[i] = rand() % (_m1-1);
		std::cout << "_m1=" << _m1 << std::endl;
		std::cout << "_a1[i]" << _a1[i] << std::endl;
		//_a2[i] = rnd.Integer<unsigned long long>(_m2-1);
		_a2[i] = rand() % (_m2-1);
	}			
}

void
CHashFunc::UHashfunc_init(
						  unsigned int t, 
						  unsigned int n, 
						  unsigned int c,
						  int32 newseed)
{	
	// Init member variables
	_t = t;
	_n = n;
	_c = c;
	
	// Get the prime number which is larger than t*n 
	unsigned long long top = _t*_n;
	
	unsigned long long p = 0;
	unsigned int mul = 1;
	do {		
		unsigned from = 100 * mul;
		p = GetPrime(top, from);
		++mul;
	} while (p == 0);
	
	_m1 = p;   	// t*n ~~ m1
	
	unsigned long long top2 = _c*_t*_n;
	unsigned long long p2 = 0;
	mul = 1;
	do {		
		unsigned from = 100 * mul;
		p2 = GetPrime(top2, from);
		++mul;
	} while (p2 == 0);
	
	_m2 = p2; // m2 > c*t*n --> to avoid double collision ==> I just use _c*top for _m2
    if(_a1 != NULL)
        delete[] _a1;
    if(_a2 != NULL)
        delete[] _a2;
	_a1 = new unsigned long long[_n];
	_a2 = new unsigned long long[_n];
	
	// generate n random numbers between 0 and m1-1
	// for hash function1 and hash function2
	// rand() % 48       
	// random number between 0~47
	//RandomLib::Random rnd(newseed);		// r created with random seed
	//std::cout << "    Random seed set to " << rnd.SeedString() << std::endl;
	
    srand(newseed);
	for (unsigned int i=0; i<_n; ++i) {
		//_a1[i] = rnd.Integer<unsigned long long>(_m1-1);
		_a1[i] = rand() % (_m1-1);
		//_a2[i] = rnd.Integer<unsigned long long>(_m2-1);
		_a2[i] = rand() % (_m2-1);
	}			
	
	//	// Init member variables
	//	_t = t;
	//	_n = n;
	////	_r = r; // hash table increase factor
	//	_c = c;
	//		
	////	std::cout << "LONG_MAX=" << LONG_MAX << std::endl;
	////	std::cout << "ULONG_MAX=" << ULONG_MAX << std::endl;
	////	std::cout << "LLONG_MAX=" << LLONG_MAX << std::endl;
	////	std::cout << "ULLONG_MAX=" << ULLONG_MAX << std::endl;	
	//
	//	// Get the prime number which is larger than t*n of the amount of 
	//	// t*n*r. If t=100, n=100 and r=0.1 then the size of hash table 
	//	// will be 10000+1000 = 11000
	//	unsigned long long top = _t*_n + int(_t*_n*_r);
	////	std::cout << "top1=" << top << std::endl;
	//	
	////	unsigned long p = GetPrime(top);
	//	unsigned long long p = 0;
	//	unsigned mul = 1;
	//	do {		
	//		unsigned from = 100 * mul;
	//		p = GetPrime(top, from);
	////		std::cout << "p1=" << p << std::endl;
	//		++mul;
	//	} while (p == 0);
	//	
	//	assert(top != 0);
	//	assert(p != 0);
	//	
	//	if (p >= ULONG_MAX) {
	//		std::cout << "ERROR: possible overflow!\n";
	//		exit(0);
	//	}
	//	
	//	// truncation possible !!    	
	//	_m1 = p;   	// t*n ~~ m1
	////	std::cout << "_m1=" << _m1 << std::endl;
	//	
	////	unsigned long top2 = _c*_t*_n;
	//	unsigned long top2 = _c*top;
	////	std::cout << "top2=" << top2 << std::endl;
	////	unsigned long p2 = GetPrime(top);
	//	unsigned long long p2 = 0;
	//	mul = 1;
	//	do {		
	//		unsigned from = 100 * mul;
	//		p2 = GetPrime(top2, from);
	////		std::cout << "p2=" << p2 << std::endl;
	//		++mul;
	//	} while (p2 == 0);
	//	
	//	assert(top2 != 0);
	//	assert(p2 != 0);
	//	
	//	if (p2 >= ULONG_MAX) {
	//		std::cout << "ERROR: possible overflow!\n";
	//		exit(0);
	//	}
	//	
	//	_m2 = p2; // m2 > c*t*n --> to avoid double collision ==> I just use _c*top for _m2
	////	std::cout << "_m2=" << _m2 << std::endl;
	//	
	////	int32 seed = time(0); 
	//	TRandomMersenne rg(newseed);	
	//	std::cout << "seed used= " << newseed << std::endl;
	//	
	//	std::fstream debug_out2;
	//	debug_out2.open("seed.txt", std::ios::out);
	//	debug_out2 << newseed << std::endl;
	//	debug_out2.close();
	//	
	//	// Debug
	////	std::cout << "_m1=" << _m1 << std::endl;
	////	std::cout << "_m2=" << _m2 << std::endl;
	////	std::cout << "LONG_MAX=" << LONG_MAX << std::endl;
	////	std::cout << "ULONG_MAX=" << ULONG_MAX << std::endl;
	////	std::cout << "LLONG_MAX=" << LLONG_MAX << std::endl;
	////	std::cout << "ULLONG_MAX=" << ULLONG_MAX << std::endl;	
	////	std::cout << "33306*564*_c=" << 33306*564*_c  << std::endl;	
	////	std::cout << "16384*4096*_c=" << 16384*4096*1000  << std::endl;	
	////	std::cout << "33306*8192*_c=" << 33306*8192*1000  << std::endl;	
	////	unsigned long temp1 = 33306*564*1000;
	////	unsigned long long temp2 = 33306*564*1000;
	////	std::cout << "ulong=" << temp1 <<  std::endl;	
	////	std::cout << "ulong long=" << temp2 <<  std::endl;	
	//////	printf("%I64d\n", temp1);
	//////	printf("%I64d\n", temp2);
	////	
	//////	uint64_t xx = 18446744073709551615;
	////	unsigned long long temp3 = 18446744073709551615;
	////	unsigned long x=33306;
	////	unsigned long y=567;
	////	unsigned long z=8192;
	////	
	////	uint64_t res = static_cast<uint64_t> (x*y*_c);
	////	if (res > std::numeric_limits<unsigned int>::max()) std::cout << "over1\n";
	////	res = static_cast<uint64_t> (x*z*_c);
	////	if (res > std::numeric_limits<unsigned int>::max()) std::cout << "over2\n";
	//	
	//	// universal hash funciton = a * b mod m
	//	// where B=(b1, b2,...,bn) is bit-string						 
	//	_a1 = new unsigned long long[_n];
	//	_a2 = new unsigned long long[_n];
	//		
	//	// generate n random numbers between 0 and m1-1
	//	// for hash function1 and hash function2
	//	// rand() % 48       
	//	// random number between 0~47
	//	for (unsigned i=0; i<_n; ++i) {
	//		_a1[i] = rg.IRandom(0,_m1-1);
	//		_a2[i] = rg.IRandom(0,_m2-1);
	//	}	
}

// Generate a prime number right after topNum
unsigned long long 
CHashFunc::GetPrime(unsigned long long topNum, unsigned from)
{
	unsigned long long primeNum=0;
	unsigned long long candidate=0;
	
	if (topNum <= 100) 
		candidate = 2;
	else
		candidate = topNum; 
	
	while (candidate <= topNum+from) {
		unsigned long long trialDivisor = 2; 
		int prime = 1;
		
		while (trialDivisor * trialDivisor <= candidate) {
			if (candidate % trialDivisor == 0) {
				prime = 0;
				break;
			}
			trialDivisor++;
		}
		if (prime) primeNum = candidate;
		
		candidate++;
	}
	
	return primeNum;
}

