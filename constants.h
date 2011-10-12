/* GiNaCRA - GiNaC Real Algebra package
 * Copyright (C) 2010  Ulrich Loup
 *
 * This file is part of GiNaCRA.
 *
 * GiNaCRA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GiNaCRA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GiNaCRA.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef GINAC_RA_CONSTANTS_H
#define GINAC_RA_CONSTANTS_H

/**
 * Collection of all constants needed within the GiNaC Real Algebra package.
 *
 *@author Ulrich Loup
 *@since 2010-11-01
 *@version 2011-05-26
 *
 * Notation is following http://www.possibility.com/Cpp/CppCodingStandard.html.
 */

#include <numeric> // provides accumulate(..)
#include <sstream>
#include <string>
#include <ginac/ginac.h>

namespace GiNaC
{
using namespace std;
// infenitesimal elements: epsilon << delta << gamma << zeta
const symbol varEpsilon("epsilon");
const symbol varDelta("delta");
const symbol varGamma("gamma");
const symbol varZeta("zeta");

// standard variable for standard polynomial objects
const symbol X("X");

/** Lexicographic ordering on expressions by the degree in the given symbol.
 */
struct ex_is_lesseq_deg: public binary_function<ex, ex, bool>
{
	symbol variable;
	ex_is_lesseq_deg(const symbol& variable) :
			variable(variable)
	{
	}

	/** Compares two expressions by their degree in the specified variable.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const symbol &a, const symbol &b) const
	{
		return a.degree(variable) <= b.degree(variable);
	}
};

/** Lexicographic ordering on the string representations of symbols.
 */
struct symbol_is_lesseq_lex: public binary_function<symbol, symbol, bool>
{
	/** Compares two expressions lexicographically by their string representations.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const symbol &a, const symbol &b) const
	{
		stringstream sA, sB;
		sA << a;
		sB << b;
		return sA.str().compare(sB.str()) <= 0;
	}
};

/** Strict lexicographic ordering on the string representations of symbols.
 */
struct symbol_is_less_lex: public binary_function<symbol, symbol, bool>
{
	/** Compares two expressions lexicographically by their string representations.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const symbol &a, const symbol &b) const
	{
		stringstream sA, sB;
		sA << a;
		sB << b;
		return sA.str().compare(sB.str()) < 0;
	}
};

/** Prototypal ordering on expressions which are terms.
 */
struct ex_is_lesseq: public binary_function<ex, ex, bool>
{
	list<symbol> variables;
	ex_is_lesseq(const list<symbol>& variables) :
			variables(variables)
	{
	}

	/** Compares two expressions lexicographically by their string representations.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const ex &a, const ex &b) const
	{
		return true;
	}

	struct get_degree: public unary_function<symbol, int>
	{
		ex p;
		get_degree(const ex &q) :
				p(q)
		{
		}
		int operator()(const symbol &x)
		{
			return p.degree(x);
		}
	};

	inline list<int> exponents(const ex &p, const list<symbol>& variables) const
	{
		int n = 0;
		for (list<symbol>::const_iterator i = variables.begin();
				i != variables.end(); ++i)
			n++;
		list<int> exponents(n, 0);
		transform(variables.begin(), variables.end(), exponents.begin(),
				get_degree(p));
		return exponents;
	}
};

/** Graded lexicographic ordering on expressions which are terms.
 */
struct ex_is_lesseq_deggrlex: public ex_is_lesseq
{
	ex_is_lesseq_deggrlex(const list<symbol>& variables) :
			ex_is_lesseq(variables)
	{
	}

	/** Compares two expressions lexicographically by their string representations.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const ex &a, const ex &b) const
	{
		list<int> e1 = exponents(a, variables);
		list<int> e2 = exponents(b, variables);
		int deg1 = accumulate(e1.begin(), e1.end(), 0);
		int deg2 = accumulate(e2.begin(), e2.end(), 0);
		if (deg1 < deg2)
			return true;
		else if (deg1 > deg2)
			return false;
		else
			return !lexicographical_compare(e1.begin(), e1.end(), e2.begin(),
					e2.end());
	}
};

/** Lexicographic degree-reverse ordering on expressions which are terms.
 */
struct ex_is_lesseq_degrevlex: public ex_is_lesseq
{
	ex_is_lesseq_degrevlex(const list<symbol>& variables) :
			ex_is_lesseq(variables)
	{
	}

	/** Compares two expressions lexicographically by their string representations.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const ex &a, const ex &b) const
	{
		list<int> e1 = exponents(a, variables);
		list<int> e2 = exponents(b, variables);
		int deg1 = accumulate(e1.begin(), e1.end(), 0);
		int deg2 = accumulate(e2.begin(), e2.end(), 0);
		if (deg1 < deg2)
			return true;
		else if (deg1 > deg2)
			return false;
		else
			return !lexicographical_compare(e1.rbegin(), e1.rend(), e2.rbegin(),
					e2.rend());
	}
};

/**
 * Relation determining whether the given expression is a monomial under the staircase, which is determined by computing the remainders dividing by the corners of a staircase.
 *@param corners corners of the staircase of a Groebner basis (unverified input)
 */
struct ex_is_under_the_staircase: public unary_function<ex, bool>
{
	list<ex> corners;
	ex_is_under_the_staircase(const list<ex> corners) :
			corners(corners)
	{
	}

	/**
	 * Tests whether the given expression is a monomial under the staircase.
	 *@param monomial (unverified input)
	 *@return true if the given expression is a monomial under the staircase, false otherwise
	 */
	bool operator()(const ex &monomial) const
	{
		for (list<ex>::const_iterator i = corners.begin(); i != corners.end();
				++i)
				{
			ex q = 0;
			divide(*i, monomial, q); // monomial divides *i by q
			if (q == 0)
				return false;
		}
		return true;
	}
};

/** Prototypal ordering on pairs of expressions which represent terms.
 */
template<class monomialOrdering>
struct expair_is_lesseq: public binary_function<pair<ex, ex> , pair<ex, ex>
		, bool>
{
	list<symbol> variables;
	expair_is_lesseq(const list<symbol>& variables) :
			variables(variables)
	{
	}

	/** Compares two expressions lexicographically by their string representations.
	 *@param a
	 *@param b
	 *@return true in case a is less or equal to b
	 */
	bool operator()(const pair<ex, ex> &a, const pair<ex, ex> &b) const
	{
		return monomialOrdering(variables)(a.second, b.second);
	}

};

} // namespace GiNaC

#endif // GINAC_RA_CONSTANTS_H
