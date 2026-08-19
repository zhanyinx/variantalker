"""
Pages package for MAFigate application.
"""

from .home import show_home_page
from .parameter_config import show_parameter_config_page
from .data_loading import show_data_loading_page

__all__ = [
    'show_home_page',
    'show_parameter_config_page',
    'show_data_loading_page'
]
