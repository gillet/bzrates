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
"""
WSGI config for mysite project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
"""

import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mysite.settings")

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
