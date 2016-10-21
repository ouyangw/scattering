///////////////////////////////////////////////////////////////////////////////
//
// Scattering: 1D scattering calculation Library
// Copyright (C) 2016 Joseph Subotnik's group, The University of Pennsylvania
//
// This file is part of Scattering.
//
// Scattering is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Scattering is distributed in the hope that it will be usefull, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Scattering. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef UTILITY_EXCEPTION_HPP
#define UTILITY_EXCEPTION_HPP
#include <exception>
#include <string>
#include <utility>
#if __cplusplus > 201100L
#define EXCEPTION_NOEXCEPT noexcept
#else
#define EXCEPTION_NOEXCEPT throw()
#endif
namespace utility
{
    class Exception : public std::exception
    {
	public:
	    Exception(const std::string& error_string, int error_code = 1):
		error_code_(error_code),
		error_string_(error_string)
        {}
#if __cplusplus > 201100L
            Exception(std::string&& error_string, int error_code = 1) noexcept:
                error_code_(std::move(error_code)),
                error_string_(std::move(error_string))
        {}
#endif
	    virtual ~Exception() EXCEPTION_NOEXCEPT
	    {}

	    const char* what() const EXCEPTION_NOEXCEPT
	    {
		return error_string_.c_str();
	    }

	    virtual int getErrorCode() const EXCEPTION_NOEXCEPT
	    {
		return error_code_;
	    }
	private:
	    int error_code_;
	    std::string error_string_;
    };
}
#endif
