/*
    GpoSolver, Library for Global Polynomial Optimization
    Copyright (C) 2014-2015 Jan Heller, <hellej1@cmp.felk.cvut.cz>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef LOGGING_H_
#define LOGGING_H_

#include <sstream>
#include <string>

namespace GpoSolver
{

    using namespace std;

    class SimpleLogger
    {
    public:
        typedef enum
        {
            INFO,
            WARNING,
            FATAL
        } LogLevel;

        static SimpleLogger &getInstance()
        {
            static SimpleLogger instance;
            return instance;
        }

        static void setLevel(const LogLevel _vlevel)
        {
            getInstance().vlevel = _vlevel;
        }

        static bool isEnabled(const LogLevel _vlevel)
        {
            return (getInstance().vlevel <= _vlevel);
        }

#ifdef NDEBUG
        static string getHeader(const LogLevel _vlevel, const char *const time)
        {
            stringstream ss;

            ss << "";

            return ss.str();
        }
#else
        static string getHeader(const LogLevel _vlevel, const char *const time, const char *const file, const int line)
        {
            stringstream ss;

            ss << "[";
            if (_vlevel == WARNING)
                ss << "WARNING ";
            else if (_vlevel == FATAL)
                ss << "FATAL ";
            ss << time << " - " << file << ":" << line << "] ";

            return ss.str();
        }
#endif

        static void log(const string s, const LogLevel _vlevel)
        {
            if (_vlevel != SimpleLogger::FATAL)
            {
#pragma omp critical(logging)
                cout << s << endl;
            }
            else
                throw runtime_error(s);
        }

    private:
        LogLevel vlevel;

        SimpleLogger() { vlevel = WARNING; };
        SimpleLogger(SimpleLogger const &);
        void operator=(SimpleLogger const &);
    };

    class ScopedLogger
    {
    public:
        ScopedLogger(const SimpleLogger::LogLevel _vlevel, const string s) : vlevel(_vlevel)
        {
            ss << s;
        }

        stringstream &stream()
        {
            return ss;
        }

        ~ScopedLogger()
        {
            SimpleLogger::getInstance().log(ss.str(), vlevel);
        }

    private:
        SimpleLogger::LogLevel vlevel;
        stringstream ss;
    };

#ifdef NDEBUG
#define LOG_INFO SimpleLogger::getInstance().isEnabled(SimpleLogger::INFO) && \
                     ScopedLogger(SimpleLogger::INFO, SimpleLogger::getInstance().getHeader(SimpleLogger::INFO, "")).stream()
#define LOG_WARN SimpleLogger::getInstance().isEnabled(SimpleLogger::WARNING) && \
                     ScopedLogger(SimpleLogger::WARNING, SimpleLogger::getInstance().getHeader(SimpleLogger::WARNING, "")).stream()
#define LOG_FATAL SimpleLogger::getInstance().isEnabled(SimpleLogger::FATAL) && \
                      ScopedLogger(SimpleLogger::FATAL, SimpleLogger::getInstance().getHeader(SimpleLogger::FATAL, "")).stream()
#else
#define LOG_INFO SimpleLogger::getInstance().isEnabled(SimpleLogger::INFO) && \
                     ScopedLogger(SimpleLogger::INFO, SimpleLogger::getInstance().getHeader(SimpleLogger::INFO, "", __FILE__, __LINE__)).stream()
#define LOG_WARN SimpleLogger::getInstance().isEnabled(SimpleLogger::WARNING) && \
                     ScopedLogger(SimpleLogger::WARNING, SimpleLogger::getInstance().getHeader(SimpleLogger::WARNING, "", __FILE__, __LINE__)).stream()
#define LOG_FATAL SimpleLogger::getInstance().isEnabled(SimpleLogger::FATAL) && \
                      ScopedLogger(SimpleLogger::FATAL, SimpleLogger::getInstance().getHeader(SimpleLogger::FATAL, "", __FILE__, __LINE__)).stream()
#endif /* NDEBUG */

}

#endif /* LOGGING_H_ */
