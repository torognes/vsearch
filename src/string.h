/*
    Copyright (C) 2015 Torbjorn Rognes

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

static char empty_string[1] = "";

class xstring
{
  char * string;
  size_t length;
  size_t alloc;
  
 public:
  
  xstring()
    {
      length = 0;
      alloc = 0;
      string = 0;
    }
  
  ~xstring()
    {
      if (alloc > 0)
        free(string);
      alloc = 0;
      string = 0;
      length = 0;
    }

  void empty()
  {
    length = 0;
  }
  
  char * get_string()
  {
    if (length > 0)
      return string;
    else
      return empty_string;
  }

  size_t get_length()
  {
    return length;
  }

  void add_c(char c)
  {
    size_t needed = 1;
    if (length + needed + 1 > alloc)
      {
        alloc = length + needed + 1;
        string = (char*) xrealloc(string, alloc);
      }
    string[length] = c;
    length += 1;
    string[length] = 0;
  }
  
  void add_d(int d)
  {
    int needed = snprintf(0, 0, "%d", d);
    if (needed < 0)
      fatal("snprintf failed");
    if (length + needed + 1 > alloc)
      {
        alloc = length + needed + 1;
        string = (char*) xrealloc(string, alloc);
      }
    sprintf(string + length, "%d", d);
    length += needed;
  }
  
  void add_s(char * s)
  {
    size_t needed = strlen(s);
    if (length + needed + 1 > alloc)
      {
        alloc = length + needed + 1;
        string = (char*) xrealloc(string, alloc);
      }
    strcpy(string + length, s);
    length += needed;
  }
};
