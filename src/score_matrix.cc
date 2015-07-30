/*
 Copyright (C) 2015 Jakob Frielingsdorf

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

 Contact: Jakob Frielingsdorf <jfrielingsdorf@gmail.com>
 */

#include "score_matrix.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "maps.h"

const char mat_blosum45[] =
    "# Entries for the BLOSUM45 matrix at a scale of ln(2)/3.0.\n\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -2 -2  0 -1 -1 -1 -1 -5\n\
R -2  7  0 -1 -3  1  0 -2  0 -3 -2  3 -1 -2 -2 -1 -1 -2 -1 -2 -1 -3  1 -1 -5\n\
N -1  0  6  2 -2  0  0  0  1 -2 -3  0 -2 -2 -2  1  0 -4 -2 -3  5 -3  0 -1 -5\n\
D -2 -1  2  7 -3  0  2 -1  0 -4 -3  0 -3 -4 -1  0 -1 -4 -2 -3  6 -3  1 -1 -5\n\
C -1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -2 -2 -3 -1 -5\n\
Q -1  1  0  0 -3  6  2 -2  1 -2 -2  1  0 -4 -1  0 -1 -2 -1 -3  0 -2  4 -1 -5\n\
E -1  0  0  2 -3  2  6 -2  0 -3 -2  1 -2 -3  0  0 -1 -3 -2 -3  1 -3  5 -1 -5\n\
G  0 -2  0 -1 -3 -2 -2  7 -2 -4 -3 -2 -2 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -5\n\
H -2  0  1  0 -3  1  0 -2 10 -3 -2 -1  0 -2 -2 -1 -2 -3  2 -3  0 -2  0 -1 -5\n\
I -1 -3 -2 -4 -3 -2 -3 -4 -3  5  2 -3  2  0 -2 -2 -1 -2  0  3 -3  4 -3 -1 -5\n\
L -1 -2 -3 -3 -2 -2 -2 -3 -2  2  5 -3  2  1 -3 -3 -1 -2  0  1 -3  4 -2 -1 -5\n\
K -1  3  0  0 -3  1  1 -2 -1 -3 -3  5 -1 -3 -1 -1 -1 -2 -1 -2  0 -3  1 -1 -5\n\
M -1 -1 -2 -3 -2  0 -2 -2  0  2  2 -1  6  0 -2 -2 -1 -2  0  1 -2  2 -1 -1 -5\n\
F -2 -2 -2 -4 -2 -4 -3 -3 -2  0  1 -3  0  8 -3 -2 -1  1  3  0 -3  1 -3 -1 -5\n\
P -1 -2 -2 -1 -4 -1  0 -2 -2 -2 -3 -1 -2 -3  9 -1 -1 -3 -3 -3 -2 -3 -1 -1 -5\n\
S  1 -1  1  0 -1  0  0  0 -1 -2 -3 -1 -2 -2 -1  4  2 -4 -2 -1  0 -2  0 -1 -5\n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1  2  5 -3 -1  0  0 -1 -1 -1 -5\n\
W -2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2  1 -3 -4 -3 15  3 -3 -4 -2 -2 -1 -5\n\
Y -2 -1 -2 -2 -3 -1 -2 -3  2  0  0 -1  0  3 -3 -2 -1  3  8 -1 -2  0 -2 -1 -5\n\
V  0 -2 -3 -3 -1 -3 -3 -3 -3  3  1 -2  1  0 -3 -1  0 -3 -1  5 -3  2 -3 -1 -5\n\
B -1 -1  5  6 -2  0  1 -1  0 -3 -3  0 -2 -3 -2  0  0 -4 -2 -3  5 -3  1 -1 -5\n\
J -1 -3 -3 -3 -2 -2 -3 -4 -2  4  4 -3  2  1 -3 -2 -1 -2  0  2 -3  4 -2 -1 -5\n\
Z -1  1  0  1 -3  4  5 -2  0 -3 -2  1 -1 -3 -1  0 -1 -2 -2 -3  1 -2  5 -1 -5\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -5\n\
* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1\n";

const char mat_blosum50[] =
    "# Entries for the BLOSUM50 matrix at a scale of ln(2)/3.0.\n\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -2 -1 -1 -5\n\
R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1 -3  0 -1 -5\n\
N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  5 -4  0 -1 -5\n\
D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  6 -4  1 -1 -5\n\
C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -2 -3 -1 -5\n\
Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0 -3  4 -1 -5\n\
E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1 -3  5 -1 -5\n\
G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -4 -2 -1 -5\n\
H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0 -3  0 -1 -5\n\
I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4  4 -3 -1 -5\n\
L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4  4 -3 -1 -5\n\
K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0 -3  1 -1 -5\n\
M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3  2 -1 -1 -5\n\
F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4  1 -4 -1 -5\n\
P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -3 -1 -1 -5\n\
S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0 -3  0 -1 -5\n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1 -1 -1 -5\n\
W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -2 -1 -5\n\
Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -1 -2 -1 -5\n\
V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -3  2 -3 -1 -5\n\
B -2 -1  5  6 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -3  6 -4  1 -1 -5\n\
J -2 -3 -4 -4 -2 -3 -3 -4 -3  4  4 -3  2  1 -3 -3 -1 -2 -1  2 -4  4 -3 -1 -5\n\
Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  1 -3  5 -1 -5\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -5\n\
* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1\n";

const char mat_blosum62[] =
    "# Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0.\n\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -1 -4\n\
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1 -2  0 -1 -4\n\
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  4 -3  0 -1 -4\n\
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4 -3  1 -1 -4\n\
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -1 -3 -1 -4\n\
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0 -2  4 -1 -4\n\
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1 -3  4 -1 -4\n\
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -4\n\
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0 -3  0 -1 -4\n\
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3  3 -3 -1 -4\n\
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4  3 -3 -1 -4\n\
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0 -3  1 -1 -4\n\
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3  2 -1 -1 -4\n\
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3  0 -3 -1 -4\n\
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -3 -1 -1 -4\n\
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0 -2  0 -1 -4\n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -1 -4\n\
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -2 -2 -1 -4\n\
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -1 -2 -1 -4\n\
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3  2 -2 -1 -4\n\
B -2 -1  4  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4 -3  0 -1 -4\n\
J -1 -2 -3 -3 -1 -2 -3 -4 -3  3  3 -3  2  0 -3 -2 -1 -2 -1  2 -3  3 -3 -1 -4\n\
Z -1  0  0  1 -3  4  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -2 -2 -2  0 -3  4 -1 -4\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4\n\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

const char mat_blosum80[] =
    "# Entries for the BLOSUM80 matrix at a scale of ln(2)/2.0.\n\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  5 -2 -2 -2 -1 -1 -1  0 -2 -2 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -2 -1 -1 -6\n\
R -2  6 -1 -2 -4  1 -1 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -4 -3 -3 -1 -3  0 -1 -6\n\
N -2 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4  5 -4  0 -1 -6\n\
D -2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -6 -4 -4  5 -5  1 -1 -6\n\
C -1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1 -4 -2 -4 -1 -6\n\
Q -1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3  0 -3  4 -1 -6\n\
E -1 -1 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3  1 -4  5 -1 -6\n\
G  0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4 -1 -5 -3 -1 -6\n\
H -2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  2 -4 -1 -4  0 -1 -6\n\
I -2 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3 -4  3 -4 -1 -6\n\
L -2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -3 -3 -2 -2 -2  1 -4  3 -3 -1 -6\n\
K -1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3 -1 -3  1 -1 -6\n\
M -1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1 -3  2 -1 -1 -6\n\
F -3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1 -4  0 -4 -1 -6\n\
P -1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -3 -1 -3 -4  8 -1 -2 -5 -4 -3 -2 -4 -2 -1 -6\n\
S  1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2  0 -3  0 -1 -6\n\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0 -1 -1 -1 -1 -6\n\
W -3 -4 -4 -6 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -5 -4 -4 11  2 -3 -5 -3 -3 -1 -6\n\
Y -2 -3 -3 -4 -3 -2 -3 -4  2 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2 -3 -2 -3 -1 -6\n\
V  0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  4 -4  2 -3 -1 -6\n\
B -2 -1  5  5 -4  0  1 -1 -1 -4 -4 -1 -3 -4 -2  0 -1 -5 -3 -4  5 -4  0 -1 -6\n\
J -2 -3 -4 -5 -2 -3 -4 -5 -4  3  3 -3  2  0 -4 -3 -1 -3 -2  2 -4  3 -3 -1 -6\n\
Z -1  0  0  1 -4  4  5 -3  0 -4 -3  1 -1 -4 -2  0 -1 -3 -3 -3  0 -3  5 -1 -6\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -6\n\
* -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1\n";

const char mat_blosum90[] =
    "# Entries for the BLOSUM90 matrix at a scale of ln(2)/2.0.\n\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *\n\
A  5 -2 -2 -3 -1 -1 -1  0 -2 -2 -2 -1 -2 -3 -1  1  0 -4 -3 -1 -2 -2 -1 -1 -6\n\
R -2  6 -1 -3 -5  1 -1 -3  0 -4 -3  2 -2 -4 -3 -1 -2 -4 -3 -3 -2 -3  0 -1 -6\n\
N -2 -1  7  1 -4  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -5 -3 -4  5 -4 -1 -1 -6\n\
D -3 -3  1  7 -5 -1  1 -2 -2 -5 -5 -1 -4 -5 -3 -1 -2 -6 -4 -5  5 -5  1 -1 -6\n\
C -1 -5 -4 -5  9 -4 -6 -4 -5 -2 -2 -4 -2 -3 -4 -2 -2 -4 -4 -2 -4 -2 -5 -1 -6\n\
Q -1  1  0 -1 -4  7  2 -3  1 -4 -3  1  0 -4 -2 -1 -1 -3 -3 -3 -1 -3  5 -1 -6\n\
E -1 -1 -1  1 -6  2  6 -3 -1 -4 -4  0 -3 -5 -2 -1 -1 -5 -4 -3  1 -4  5 -1 -6\n\
G  0 -3 -1 -2 -4 -3 -3  6 -3 -5 -5 -2 -4 -5 -3 -1 -3 -4 -5 -5 -2 -5 -3 -1 -6\n\
H -2  0  0 -2 -5  1 -1 -3  8 -4 -4 -1 -3 -2 -3 -2 -2 -3  1 -4 -1 -4  0 -1 -6\n\
I -2 -4 -4 -5 -2 -4 -4 -5 -4  5  1 -4  1 -1 -4 -3 -1 -4 -2  3 -5  3 -4 -1 -6\n\
L -2 -3 -4 -5 -2 -3 -4 -5 -4  1  5 -3  2  0 -4 -3 -2 -3 -2  0 -5  4 -4 -1 -6\n\
K -1  2  0 -1 -4  1  0 -2 -1 -4 -3  6 -2 -4 -2 -1 -1 -5 -3 -3 -1 -3  1 -1 -6\n\
M -2 -2 -3 -4 -2  0 -3 -4 -3  1  2 -2  7 -1 -3 -2 -1 -2 -2  0 -4  2 -2 -1 -6\n\
F -3 -4 -4 -5 -3 -4 -5 -5 -2 -1  0 -4 -1  7 -4 -3 -3  0  3 -2 -4  0 -4 -1 -6\n\
P -1 -3 -3 -3 -4 -2 -2 -3 -3 -4 -4 -2 -3 -4  8 -2 -2 -5 -4 -3 -3 -4 -2 -1 -6\n\
S  1 -1  0 -1 -2 -1 -1 -1 -2 -3 -3 -1 -2 -3 -2  5  1 -4 -3 -2  0 -3 -1 -1 -6\n\
T  0 -2  0 -2 -2 -1 -1 -3 -2 -1 -2 -1 -1 -3 -2  1  6 -4 -2 -1 -1 -2 -1 -1 -6\n\
W -4 -4 -5 -6 -4 -3 -5 -4 -3 -4 -3 -5 -2  0 -5 -4 -4 11  2 -3 -6 -3 -4 -1 -6\n\
Y -3 -3 -3 -4 -4 -3 -4 -5  1 -2 -2 -3 -2  3 -4 -3 -2  2  8 -3 -4 -2 -3 -1 -6\n\
V -1 -3 -4 -5 -2 -3 -3 -5 -4  3  0 -3  0 -2 -3 -2 -1 -3 -3  5 -4  1 -3 -1 -6\n\
B -2 -2  5  5 -4 -1  1 -2 -1 -5 -5 -1 -4 -4 -3  0 -1 -6 -4 -4  5 -5  0 -1 -6\n\
J -2 -3 -4 -5 -2 -3 -4 -5 -4  3  4 -3  2  0 -4 -3 -2 -3 -2  1 -5  4 -4 -1 -6\n\
Z -1  0 -1  1 -5  5  5 -3  0 -4 -4  1 -2 -4 -2 -1 -1 -4 -3 -3  0 -4  5 -1 -6\n\
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -6\n\
* -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1\n";

const char mat_pam30[] =
    "#\n\
# This matrix was produced by \"pam\" Version 1.0.6 [28-Jul-93]\n\
#\n\
# PAM 30 substitution matrix, scale = ln(2)/2 = 0.346574\n\
#\n\
# Expected score = -5.06, Entropy = 2.57 bits\n\
#\n\
# Lowest score = -17, Highest score = 13\n\
#\n\
    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *\n\
A   6  -7  -4  -3  -6  -4  -2  -2  -7  -5  -6  -7  -5  -8  -2   0  -1 -13  -8  -2  -3  -3  -3 -17\n\
R  -7   8  -6 -10  -8  -2  -9  -9  -2  -5  -8   0  -4  -9  -4  -3  -6  -2 -10  -8  -7  -4  -6 -17\n\
N  -4  -6   8   2 -11  -3  -2  -3   0  -5  -7  -1  -9  -9  -6   0  -2  -8  -4  -8   6  -3  -3 -17\n\
D  -3 -10   2   8 -14  -2   2  -3  -4  -7 -12  -4 -11 -15  -8  -4  -5 -15 -11  -8   6   1  -5 -17\n\
C  -6  -8 -11 -14  10 -14 -14  -9  -7  -6 -15 -14 -13 -13  -8  -3  -8 -15  -4  -6 -12 -14  -9 -17\n\
Q  -4  -2  -3  -2 -14   8   1  -7   1  -8  -5  -3  -4 -13  -3  -5  -5 -13 -12  -7  -3   6  -5 -17\n\
E  -2  -9  -2   2 -14   1   8  -4  -5  -5  -9  -4  -7 -14  -5  -4  -6 -17  -8  -6   1   6  -5 -17\n\
G  -2  -9  -3  -3  -9  -7  -4   6  -9 -11 -10  -7  -8  -9  -6  -2  -6 -15 -14  -5  -3  -5  -5 -17\n\
H  -7  -2   0  -4  -7   1  -5  -9   9  -9  -6  -6 -10  -6  -4  -6  -7  -7  -3  -6  -1  -1  -5 -17\n\
I  -5  -5  -5  -7  -6  -8  -5 -11  -9   8  -1  -6  -1  -2  -8  -7  -2 -14  -6   2  -6  -6  -5 -17\n\
L  -6  -8  -7 -12 -15  -5  -9 -10  -6  -1   7  -8   1  -3  -7  -8  -7  -6  -7  -2  -9  -7  -6 -17\n\
K  -7   0  -1  -4 -14  -3  -4  -7  -6  -6  -8   7  -2 -14  -6  -4  -3 -12  -9  -9  -2  -4  -5 -17\n\
M  -5  -4  -9 -11 -13  -4  -7  -8 -10  -1   1  -2  11  -4  -8  -5  -4 -13 -11  -1 -10  -5  -5 -17\n\
F  -8  -9  -9 -15 -13 -13 -14  -9  -6  -2  -3 -14  -4   9 -10  -6  -9  -4   2  -8 -10 -13  -8 -17\n\
P  -2  -4  -6  -8  -8  -3  -5  -6  -4  -8  -7  -6  -8 -10   8  -2  -4 -14 -13  -6  -7  -4  -5 -17\n\
S   0  -3   0  -4  -3  -5  -4  -2  -6  -7  -8  -4  -5  -6  -2   6   0  -5  -7  -6  -1  -5  -3 -17\n\
T  -1  -6  -2  -5  -8  -5  -6  -6  -7  -2  -7  -3  -4  -9  -4   0   7 -13  -6  -3  -3  -6  -4 -17\n\
W -13  -2  -8 -15 -15 -13 -17 -15  -7 -14  -6 -12 -13  -4 -14  -5 -13  13  -5 -15 -10 -14 -11 -17\n\
Y  -8 -10  -4 -11  -4 -12  -8 -14  -3  -6  -7  -9 -11   2 -13  -7  -6  -5  10  -7  -6  -9  -7 -17\n\
V  -2  -8  -8  -8  -6  -7  -6  -5  -6   2  -2  -9  -1  -8  -6  -6  -3 -15  -7   7  -8  -6  -5 -17\n\
B  -3  -7   6   6 -12  -3   1  -3  -1  -6  -9  -2 -10 -10  -7  -1  -3 -10  -6  -8   6   0  -5 -17\n\
Z  -3  -4  -3   1 -14   6   6  -5  -1  -6  -7  -4  -5 -13  -4  -5  -6 -14  -9  -6   0   6  -5 -17\n\
X  -3  -6  -3  -5  -9  -5  -5  -5  -5  -5  -6  -5  -5  -8  -5  -3  -4 -11  -7  -5  -5  -5  -5 -17\n\
* -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17 -17   1\n";

const char mat_pam70[] =
    "#\n\
# This matrix was produced by \"pam\" Version 1.0.6 [28-Jul-93]\n\
#\n\
# PAM 70 substitution matrix, scale = ln(2)/2 = 0.346574\n\
#\n\
# Expected score = -2.77, Entropy = 1.60 bits\n\
#\n\
# Lowest score = -11, Highest score = 13\n\
#\n\
    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *\n\
A   5  -4  -2  -1  -4  -2  -1   0  -4  -2  -4  -4  -3  -6   0   1   1  -9  -5  -1  -1  -1  -2 -11\n\
R  -4   8  -3  -6  -5   0  -5  -6   0  -3  -6   2  -2  -7  -2  -1  -4   0  -7  -5  -4  -2  -3 -11\n\
N  -2  -3   6   3  -7  -1   0  -1   1  -3  -5   0  -5  -6  -3   1   0  -6  -3  -5   5  -1  -2 -11\n\
D  -1  -6   3   6  -9   0   3  -1  -1  -5  -8  -2  -7 -10  -4  -1  -2 -10  -7  -5   5   2  -3 -11\n\
C  -4  -5  -7  -9   9  -9  -9  -6  -5  -4 -10  -9  -9  -8  -5  -1  -5 -11  -2  -4  -8  -9  -6 -11\n\
Q  -2   0  -1   0  -9   7   2  -4   2  -5  -3  -1  -2  -9  -1  -3  -3  -8  -8  -4  -1   5  -2 -11\n\
E  -1  -5   0   3  -9   2   6  -2  -2  -4  -6  -2  -4  -9  -3  -2  -3 -11  -6  -4   2   5  -3 -11\n\
G   0  -6  -1  -1  -6  -4  -2   6  -6  -6  -7  -5  -6  -7  -3   0  -3 -10  -9  -3  -1  -3  -3 -11\n\
H  -4   0   1  -1  -5   2  -2  -6   8  -6  -4  -3  -6  -4  -2  -3  -4  -5  -1  -4   0   1  -3 -11\n\
I  -2  -3  -3  -5  -4  -5  -4  -6  -6   7   1  -4   1   0  -5  -4  -1  -9  -4   3  -4  -4  -3 -11\n\
L  -4  -6  -5  -8 -10  -3  -6  -7  -4   1   6  -5   2  -1  -5  -6  -4  -4  -4   0  -6  -4  -4 -11\n\
K  -4   2   0  -2  -9  -1  -2  -5  -3  -4  -5   6   0  -9  -4  -2  -1  -7  -7  -6  -1  -2  -3 -11\n\
M  -3  -2  -5  -7  -9  -2  -4  -6  -6   1   2   0  10  -2  -5  -3  -2  -8  -7   0  -6  -3  -3 -11\n\
F  -6  -7  -6 -10  -8  -9  -9  -7  -4   0  -1  -9  -2   8  -7  -4  -6  -2   4  -5  -7  -9  -5 -11\n\
P   0  -2  -3  -4  -5  -1  -3  -3  -2  -5  -5  -4  -5  -7   7   0  -2  -9  -9  -3  -4  -2  -3 -11\n\
S   1  -1   1  -1  -1  -3  -2   0  -3  -4  -6  -2  -3  -4   0   5   2  -3  -5  -3   0  -2  -1 -11\n\
T   1  -4   0  -2  -5  -3  -3  -3  -4  -1  -4  -1  -2  -6  -2   2   6  -8  -4  -1  -1  -3  -2 -11\n\
W  -9   0  -6 -10 -11  -8 -11 -10  -5  -9  -4  -7  -8  -2  -9  -3  -8  13  -3 -10  -7 -10  -7 -11\n\
Y  -5  -7  -3  -7  -2  -8  -6  -9  -1  -4  -4  -7  -7   4  -9  -5  -4  -3   9  -5  -4  -7  -5 -11\n\
V  -1  -5  -5  -5  -4  -4  -4  -3  -4   3   0  -6   0  -5  -3  -3  -1 -10  -5   6  -5  -4  -2 -11\n\
B  -1  -4   5   5  -8  -1   2  -1   0  -4  -6  -1  -6  -7  -4   0  -1  -7  -4  -5   5   1  -2 -11\n\
Z  -1  -2  -1   2  -9   5   5  -3   1  -4  -4  -2  -3  -9  -2  -2  -3 -10  -7  -4   1   5  -3 -11\n\
X  -2  -3  -2  -3  -6  -2  -3  -3  -3  -3  -4  -3  -3  -5  -3  -1  -2  -7  -5  -2  -2  -3  -3 -11\n\
* -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11 -11   1\n";

const char mat_pam250[] =
    "#\n\
# This matrix was produced by \"pam\" Version 1.0.6 [28-Jul-93]\n\
#\n\
# PAM 250 substitution matrix, scale = ln(2)/3 = 0.231049\n\
#\n\
# Expected score = -0.844, Entropy = 0.354 bits\n\
#\n\
# Lowest score = -8, Highest score = 17\n\
#\n\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X\n\
A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0  0  0  0\n\
R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2 -1  0 -1\n\
N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2  2  1  0\n\
D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2  3  3 -1\n\
C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2 -4 -5 -3\n\
Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2  1  3 -1\n\
E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2  3  3 -1\n\
G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1  0  0 -1\n\
H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2  1  2 -1\n\
I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4 -2 -2 -1\n\
L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2 -3 -3 -1\n\
K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2  1  0 -1\n\
M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2 -2 -2 -1\n\
F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1 -4 -5 -2\n\
P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1 -1  0 -1\n\
S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1  0  0  0\n\
T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0  0 -1  0\n\
W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6 -5 -6 -4\n\
Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2 -3 -4 -2\n\
V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4 -2 -2 -1\n\
B  0 -1  2  3 -4  1  3  0  1 -2 -3  1 -2 -4 -1  0  0 -5 -3 -2  3  2 -1\n\
Z  0  0  1  3 -5  3  3  0  2 -2 -3  0 -2 -5  0  0 -1 -6 -4 -2  2  3 -1\n\
X  0 -1  0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1  0  0 -4 -2 -1 -1 -1 -1\n";

ScoreMatrix ScoreMatrix::instance = ScoreMatrix();

void ScoreMatrix::init(const char* amino_acid_matrix)
{
  constant_scoring = 0;
  mode = MATRIX_MODE_AA;

  prepare_matrices(MATRIX_MODE_AA);

  mat_init_buildin(amino_acid_matrix);

  finalize_matrices();
}

void ScoreMatrix::init(const int match, const int mismatch, int sequence_mode)
{
  constant_scoring = 1;
  mode = sequence_mode;

  prepare_matrices(sequence_mode);

  mat_init_constant_scoring(match, mismatch);

  finalize_matrices();
}

ScoreMatrix::~ScoreMatrix()
{
  mat_free();
}

/**
 * Allocates memory for storing the scoring matrices.
 *
 * score_matrix_64 is initialised with -1.
 */
void ScoreMatrix::prepare_matrices(int sequence_mode)
{
  if (sequence_mode==MATRIX_MODE_NUC)
    dimension = 16;
  else
    dimension = 32;

  // free first, before we allocate new memory
  mat_free();

  score_matrix_16 = (CELL*)xmalloc(dimension*dimension*sizeof(CELL));
  score_matrix_64 = (long*)xmalloc(dimension*dimension*sizeof(long));

  memset(score_matrix_64, -1, dimension*dimension*sizeof(long));
}

/**
 * Copies the values read into score_matrix_64 into the other representations
 * of the scoring matrices.
 */
void ScoreMatrix::finalize_matrices()
{
  int a, b;
  long sc, lo, hi;

  hi = -100;
  lo = 100;

  for (a = 0; a<dimension; a++)
    {
      for (b = 0; b<dimension; b++)
        {
          sc = get64(a, b);
          if (sc<lo)
            lo = sc;
          if (sc>hi)
            hi = sc;
        }
    }

  for (a = 0; a<dimension; a++)
    {
      for (b = 0; b<dimension; b++)
        {
          set16(a, b, get64(a, b));
        }
    }
}

/**
 * Reads a line into score_matrix_64
 *
 * @param line      the line to read in
 * @param symbols   nr of read symbols
 * @param order     TODO ???
 */
void ScoreMatrix::read_line(char line[LINE_MAX], int* symbols, char* order)
{
  char *p, *q;

  p = line;
  char c = *p++;

  switch (c)
    {
    case '\n':
      case '#':
      /* ignore blank lines and comments starting with # */

      break;
    case '\t':
      case ' ':
      /* read order of symbols, copy non-whitespace chars */

      q = order;

      while ((c = *p++))
        {
          if (strchr(" \t\n", c)==NULL)
            {
              *q++ = chrmap_aa_5bit[(int)c];
              (*symbols)++;
            }
        }

      break;
    default:
      {
        /* ordinary lines */
        int read;

        int a = chrmap_aa_5bit[(int)c];
        for (int i = 0; i<*symbols; i++)
          {
            long sc;

            if (sscanf(p, "%ld%n", &sc, &read)==0)
              {
                fatal("Problem parsing score matrix file.");
              }

            char b = order[i];

            if ((a>=0)&&(b>=0)&&(a<dimension)&&(b<dimension))
              {
                set64(a, b, sc);
              }

            p += read;
          }
        break;
      }
    }
}

/**
 * Releases the memory allocated for the scoring matrices.
 */
void ScoreMatrix::mat_free()
{
  if (score_matrix_16)
    {
      free(score_matrix_16);
      score_matrix_16 = NULL;
    }

  if (score_matrix_64)
    {
      free(score_matrix_64);
      score_matrix_64 = NULL;
    }
}

void ScoreMatrix::mat_init_from_string(const char * matrix)
{
  char line[LINE_MAX];
  char order[LINE_MAX];

  int symbols;

  char * s = (char*)matrix;

  if (!s)
    fatal("Cannot read score matrix string.");

  symbols = 0;

  while (*s)
    {
      char * nextline = strchr(s, '\n');
      size_t linelen;
      if (nextline)
        linelen = nextline-s;
      else
        linelen = strlen(s);

      strncpy(line, s, linelen);
      line[linelen] = 0;

      read_line(line, &symbols, order);

      if (nextline)
        s = nextline+1;
      else
        s = s+strlen(s);
    }
}

void ScoreMatrix::mat_init_from_file(const char * matrix)
{
  // reading of the file
  char line[LINE_MAX];
  char order[LINE_MAX];

  int symbols;

  FILE * fp = fopen(matrix, "r");

  if (!fp)
    fatal("Cannot open score matrix file.");

  symbols = 0;

  while (fgets(line, LINE_MAX, fp)!=NULL)
    {
      read_line(line, &symbols, order);
    }

  fclose(fp);
}

void ScoreMatrix::mat_init_constant_scoring(const int match, const int mismatch)
{
  for (int a = 1; a<dimension; a++)
    {
      for (int b = 1; b<dimension; b++)
        {
          set64(a, b, ((a==b) ? match : mismatch));
        }
    }
}

void ScoreMatrix::mat_init_buildin(const char* matrixname)
{
  if (strcasecmp(matrixname, BLOSUM45)==0)
    mat_init_from_string(mat_blosum45);
  else if (strcasecmp(matrixname, BLOSUM50)==0)
    mat_init_from_string(mat_blosum50);
  else if (strcasecmp(matrixname, BLOSUM62)==0)
    mat_init_from_string(mat_blosum62);
  else if (strcasecmp(matrixname, BLOSUM80)==0)
    mat_init_from_string(mat_blosum80);
  else if (strcasecmp(matrixname, BLOSUM90)==0)
    mat_init_from_string(mat_blosum90);
  else if (strcasecmp(matrixname, PAM30)==0)
    mat_init_from_string(mat_pam30);
  else if (strcasecmp(matrixname, PAM70)==0)
    mat_init_from_string(mat_pam70);
  else if (strcasecmp(matrixname, PAM250)==0)
    mat_init_from_string(mat_pam250);
  else
    fatal("Unknown matrix: %s", matrixname);
}

void ScoreMatrix::dump_matrix()
{
  printf("     ");
  for (int i = 0; i<dimension; i++)
    printf("%2d", i);
  printf("\n");
  printf("     ");
  for (int i = 0; i<dimension; i++)
    printf(" %c", sym_aa_5bit[i]);
  printf("\n");
  for (int i = 0; i<dimension; i++)
    {
      printf("%2d %c ", i, sym_aa_5bit[i]);
      for (int j = 0; j<dimension; j++)
        {
          printf("%2ld", get64(i, j));
        }
      printf("\n");
    }
}
