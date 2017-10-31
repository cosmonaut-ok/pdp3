## Copyright (C) 2011 Philip
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## getframe

## Author: Philip <address@hidden>
## Created: 2011-12-14

function [ ret ] = getframe (h)

  print (h, "/tmp/tmp.fig", "-dppm");
  ret = im2double (imread ("/tmp/tmp.fig"));
  ## Truncate to even size to accomodate addframe()
  if (mod (size (ret, 1), 2) > 0); ret = ret(2:end, :, :); endif
  if (mod (size (ret, 2), 2) > 0); ret = ret(:, 2:end, :); endif

endfunction

