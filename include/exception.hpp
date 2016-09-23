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
