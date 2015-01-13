
import os, sys
sys.path.append('/var/www/bzrates')
sys.path.append('/var/www/bzrates/mysite')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mysite.settings")

from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
