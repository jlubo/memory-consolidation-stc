/********************************************
 *** Useful generic functions and classes ***
 ********************************************/

/*** Copyright 2017-2022 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

// namespace std is required

#include <vector> // for this header, add -std=c++0x to the command line
#include <chrono> // for dateStr() and getClockSeed()
#include <limits> // for maximum/minimum magnitude of integers
#include <cmath>
#include <cstring>
#include <sys/ioctl.h> // for getLinebreak()
#include <errno.h>
#include <memory> // for unique pointers


/*** concat ***
 * Concatenates two character arrays/strings and returns a string *
 * - str1: the first character array
 * - str2: the second character array
 * - return: the concatenated character array in a new string variable */
string concat(const char *str1, const char *str2)
{
	string s = string(str1) + string(str2);
	return s;
}
string concat(string str1, const char *str2)
{
	string s = str1 + string(str2);
	return s;
}
string concat(const char *str1, string str2)
{
	string s = string(str1) + str2;
	return s;
}


/*** pow2 ***
 * Returns the given number to the power 2 *
 * - number: the number of which the power shall be determined
 * - return: the number by the power 2 */
inline double pow2 (double number)
{
	return number*number;
}
inline int pow2 (int number)
{
	return number*number;
}

/*** pow_int ***
 * Returns the given number to the given power, all-integer function *
 * - number: the number of which the power shall be determined
 * - p: the power
 * - return: the number by the power */
inline int pow_int (int number, int p)
{
	int ret = 1;

	for (int i=0; i<p; i++)
		ret *= number;

	return ret;
}

/*** dtos ***
 * Converts a double variable to a string, rounded to the given number of internal decimal places, *
 * using the C method sprintf *
 * - num: the double floating point number to be converted *
 * - n [optional]: the number of internal decimal places (has to be greater than or equal to 0) *
 * - dont_allow_zero [optional]: specifies if zeros shall be replaced by "nan" values *
 * - return: the string containing the rounded floating point number */
string dtos(double num, int n = 0, bool dont_allow_zero = false)
{
	string s("nan"); // the string to be returned
	char *buf; // the buffer the converted number is written to
	char *format; // the format string for sprintf
	int prepoint; // number of places before point, has to be at least 1
	int dn = (n > 0 ? (int) floor(log10(n)) + 1 : 1); // decimal places of the variable n (needed to generate the format string)

	if (!std::isnan(num))
	{
		if (num <= 0.0 || floor(log10(num)) < 0.0)
			prepoint = 1;
		else
			prepoint = (int) floor(log10(num)) + 1;

		buf = new char[prepoint + n + 2]; // decimal places plus internal decimal places of num plus 2 for point and zero-terminator
		format = new char[dn + 4]; // decimal places of n plus 4 for the characters "%.f\0"

		sprintf(format, "%%.%df", n); // create format string
		sprintf(buf, (const char*) format, num); // create character array from double
	
		if (!(num == 0. && dont_allow_zero)) // if value will NOT be "nan"
			s.assign(buf); // create buffer

		delete[] buf;
		delete[] format;
	}

	return s;

	/* alternative (does not strip trailing zeros, though)

	if (!(std::isnan(num) || (num == 0. && dont_allow_zero))) // if value will NOT be "nan"
	{
		int factor = pow_int(10, n);
		int num_as_int = num * factor;

		return to_string(double(num_as_int) / factor);
	}

	return string("nan");*/

}


/*** dateStr ***
 * Returns a string containing a time stamp of 18 char's, *
 * optionally added to the start of another string *
 * - str [optional]: a character array to which the time stamp will be added  *
 * - fixdate [optional]: true if this time shall be fixed, false by default *
 * - return: the date string plus the specified string, if stated */
string dateStr(string str = "", bool fixdate = false) 
{
	const int len = 19; // time stamp length: 18 characters (plus terminating \0)
	static string datestr("");

	if (fixdate == true) // if this date is to be fixed, update date
	{
		char buf[len];
		time_t t = time(NULL);
		struct tm *now  = localtime(&t);
		strftime(buf, len, "%y-%m-%d_%H-%M-%S", now); 
		datestr.assign(buf);
	}

	return datestr + str;
}

/*** copyFile ***
 * Copies the content of a file from a source file to a target file *
 * - src: file name of the source file
 * - dest: file name of the target file
 * - return: false if one of the files could not be opened, true if succesful */
bool copyFile(string src, string dest)
{
	ifstream ssrc;
	ofstream sdest;
	
	ssrc.open(src, ios::binary);
	if (!ssrc.is_open())
		return false;

	sdest.open(dest, ios::binary);
	if (!sdest.is_open())
		return false;
	
	sdest << ssrc.rdbuf();

	ssrc.close();
	sdest.close();
	
	return true;
}

/*** extractFileFromBinary ***
 * Extracts a file that was attached to the binary by the linker *
 * - filename: name of the file to be extracted *
 * - start_pos: start position for reading *
 * - end_pos: end position for reading */
void extractFileFromBinary(string filename, uint8_t* start_pos, uint8_t* end_pos)
{
	int i = 0;
	ofstream cf(filename, ofstream::binary);

	if (!cf.is_open())
		throw runtime_error(string("Attached file ") + string(filename) + string(" could not be extracted."));

	while (&start_pos[i] != &end_pos[0])
	{
		cf.put(start_pos[i++]);
	}

	cf.close();
}

/*** writeRunFile ***
 * Writes a file containing the arguments that were passed to the program via commandline *
 * - filename: name of the script file to be created *
 * - argc: the number of arguments *
 * - argv: array containing the arguments */
void writeRunFile(string filename, int argc, char** argv)
{
	ofstream runf(filename);

	if (!runf.is_open())
		throw runtime_error(string("\'Run\' file could not be created."));

	for(int i=0; i<argc; i++)
	{
		if (strstr(argv[i], " ") != NULL) // add quotation marks if a space occurs in argv[i]
		{
			char *arg = new char [strlen(argv[i])+1];
			strcpy(arg, argv[i]);


			char *pt = strstr(arg, "=");
			if (pt != NULL)
			{
				pt[0] = '\0';
				pt++;
				runf << arg << "=\"" << pt << "\" ";
			}
			else
			{
				runf << "\"" << arg << "\" ";
			}

			delete[] arg;
			
		}
		else
			runf << argv[i] << " ";
	}

	runf.close();
}

/*** stod_acc ***
 * String to double conversion alternative to stod *
 * - str: a string containing only a floating point number, must contain a point
 * - return: the string converted into a double */
double stod_acc(const string &str) 
{
	int place = str.find(".")-1; // decimal power
	double d = 0.0;
	
	for (int i=0; i<str.length(); i++)
	{
		double base = pow(10, place);

		if (str[i] == '1')
			d += 1.*base;
		else if (str[i] == '2')
			d += 2.*base;
		else if (str[i] == '3')
			d += 3.*base;
		else if (str[i] == '4')
			d += 4.*base;
		else if (str[i] == '5')
			d += 5.*base;
		else if (str[i] == '6')
			d += 6.*base;
		else if (str[i] == '7')
			d += 7.*base;
		else if (str[i] == '8')
			d += 8.*base;
		else if (str[i] == '9')
			d += 9.*base;
		else if (str[i] == '.')
			continue;
		place--;
	}
	return d;
}

/*** min ***
 * Finds the minimum value of a double array/a double vector/two double values *
 * - obj:  * - obj: the data object mentioned above
 * - len: the number of array elements
 * - return: the minimum value */
double min (double *obj, int len)
{
	double min_n = (len > 0) ? obj[0] : std::numeric_limits<double>::lowest();
	for (int i=1; i < len; i++)
	{
		if (obj[i] < min_n)
			min_n = obj[i];
	}
	return min_n;
}
double min (double a, double b)
{
	return (a < b ? a : b);
}

/*** max ***
 * Finds the maximum value of a double array/a double vector/two double values *
 * - obj: the data object mentioned above
 * - len: the number of array elements
 * - return: the maximum value */
double max (double *obj, int len)
{
	double max_n = (len > 0) ? obj[0] : std::numeric_limits<double>::max();
	for (int i=1; i < len; i++)
	{
		if (obj[i] > max_n)
			max_n = obj[i];
	}
	return max_n;
}
double max (vector<double> obj)
{
	double max_n = (obj.size() > 0) ? obj.front() : std::numeric_limits<double>::max();
	for (int i=1; i < obj.size(); i++)
	{
		if (obj.at(i) > max_n)
		{
			max_n = obj.at(i);
		}
	}
	return max_n;
}
double max (double a, double b)
{
	return (a > b ? a : b);
}

/*** sgn ***
 * Return the sign of a double *
 * - number: the number of which the sign shall be determined
 * - return: an integer -1 or +1 */
inline signed int sgn (double number)
{
	if (number >= 0.)
		return +1;
	else
		return -1;
}


/*** isCloser ***
 * Determines if a is closer to b than c is to b *
 * - a: a floating point number
 * - b: a floating point number
 * - return: true, if a is closer to b, false otherwise */
bool isCloser(double a, double b, double c)
{
	if ( abs(a - b) <= abs(b - c) )
		return true;
	else
		return false;
}

/*** timeMeasure ***
 * Starts or stops a time measurement (accuracy is one second) *
 * - start: boolean which is true when intending to start and false when intending to stop measurement
 * - return: the seconds elapsed */
int timeMeasure (bool start)
{
	static int start_time;
	if (start == true) 
	{
		start_time = time(NULL);
		return 0;
	}
	else
	{
		return (time(NULL) - start_time);
	}
}

/*** getMean ***
 * Calculates the mean value of a part of an array of double values *
 * - array: array of double values *
 * - start: starting point in array *
 * - len: number of array elements to be considered *
 * - return: the mean of all considered numbers */
double getMean(const double *array, const int start, const int len)
{
	double m = 0.0;
	for(int i=start; i<start+len; i++)
	{
		m += array[i];
	}
	return (m / double(len));
}

/*** getSeparator ***
 * Determines the length of a line in console/terminal and creates a line separator *
 * - return: the separator string */
string getSeparator()
{
	string sep_char("-"); // character from which the separator is created
	string separator("");
	
	struct winsize w;
	ioctl(0, TIOCGWINSZ, &w);
	
	for (int i=0; i<w.ws_col; i++)
		separator += sep_char;

	return separator;
}

/*** showChDirErrMessage ***
 * Prints the last errno message out of the possible chdir() errors via cout */
void showChDirErrMessage()
{
	cout << "Error while changing directory: ";
	switch(errno)
	{
		case EACCES:
			cout << "EACCES" << endl;
			break;
		case ENAMETOOLONG:
			cout << "ENAMETOOLONG" << endl;
			break;
		case ENOENT:
			cout << "ENOENT" << endl;
			break;
		case ENOTDIR:
			cout << "ENOTDIR" << endl;
			break;
		case ELOOP:
			cout << "ELOOP" << endl;
			break;
	}
}

/*** getClockSeed ***
 * Gets a random generator seed from the computer clock, guarantees not to return *
 * the same seed in two subsequent calls (very important!) *
 * - return: clock counts since epoch of the clock */
static unsigned int getClockSeed()
{
	static int last_seed;
	int seed;

	while ( (seed = chrono::system_clock::now().time_since_epoch().count()) == last_seed ) {}
	last_seed = seed;

	return seed;
}

/*** step ***
 * Heaviside step function, returns 1 for x >= 0 and 0 else *
 * (the alternative step_alt returns 1 for x > 0 and 0 else) *
 * - double x: argument for the step function 
 * - return: value of the Heaviside step function */
inline int step(double x)
{
	if (x >= 0.)
		return 1;
	else
		return 0;
}

/*** step_alt ***
 * Heaviside step function, returns 1 for x > 0 and 0 else *
 * - double x: argument for the step function 
 * - return: value of the Heaviside step function */
inline int step_alt(double x)
{
	if (x > 0.)
		return 1;
	else
		return 0;
}


