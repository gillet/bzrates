#!/usr/bin/env python
"""
==================================================================================
bz-rates: a web-tool to accurately estimate mutation rates from fluctuation assays
==================================================================================
Copyright (C) 2015 Alexandre Gillet-Markowska

This file is part of bz-rates.

bz-rates is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

bz-rates is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with bz-rates.  If not, see <http://www.gnu.org/licenses/>.
==================================================================================
Contact: Alexandre Gillet-Markowska - alexandre(dot)gillet(at)yahoo(dot)fr
==================================================================================
"""
import os
import sys

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mysite.settings")

    from django.core.management import execute_from_command_line

    execute_from_command_line(sys.argv)
