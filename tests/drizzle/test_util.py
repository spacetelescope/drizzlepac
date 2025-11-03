import os
import pytest
from unittest.mock import patch
from drizzlepac import util


class TestGetEnvvarSwitch:
    """Test suite for the get_envvar_switch function."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.test_envvar = 'TEST_ENVVAR_SWITCH'
        # Clean up any existing test environment variable
        if self.test_envvar in os.environ:
            del os.environ[self.test_envvar]
    
    def teardown_method(self):
        """Clean up after each test."""
        # Clean up test environment variable
        if self.test_envvar in os.environ:
            del os.environ[self.test_envvar]
    
    def test_envvar_not_set_returns_default_true(self, capsys):
        """Test that function returns default value when environment variable is not set."""
        result = util.get_envvar_switch(self.test_envvar, True, 'test setting')
        
        assert result is True
        captured = capsys.readouterr()
        assert f"ENVVAR {self.test_envvar} not found, setting test setting to default of True." in captured.out
    
    def test_envvar_not_set_returns_default_false(self, capsys):
        """Test that function returns default value when environment variable is not set."""
        result = util.get_envvar_switch(self.test_envvar, False, 'test setting')
        
        assert result is False
        captured = capsys.readouterr()
        assert f"ENVVAR {self.test_envvar} not found, setting test setting to default of False." in captured.out
    
    def test_envvar_not_set_empty_description(self, capsys):
        """Test with empty description."""
        result = util.get_envvar_switch(self.test_envvar, False, '')
        
        assert result is False
        captured = capsys.readouterr()
        assert f"ENVVAR {self.test_envvar} not found, setting to default of False." in captured.out
    
    def test_envvar_not_set_no_description(self, capsys):
        """Test with no description parameter."""
        result = util.get_envvar_switch(self.test_envvar, False)
        
        assert result is False
        captured = capsys.readouterr()
        assert f"ENVVAR {self.test_envvar} not found, setting to default of False." in captured.out
    
    @pytest.mark.parametrize("env_value,expected", [
        ("true", True),
        ("True", True),
        ("TRUE", True),
        ("yes", True),
        ("YES", True),
        ("on", True),
        ("ON", True),
        ("false", False),
        ("False", False),
        ("FALSE", False),
        ("no", False),
        ("NO", False),
        ("off", False),
        ("OFF", False),
    ])
    def test_valid_envvar_values(self, env_value, expected, capsys):
        """Test all valid environment variable values."""
        os.environ[self.test_envvar] = env_value
        
        result = util.get_envvar_switch(self.test_envvar, False, 'test setting')
        
        assert result is expected
        captured = capsys.readouterr()
        assert f"ENVVAR {self.test_envvar} found, setting test setting to {expected}." in captured.out
    
    def test_envvar_with_whitespace(self, capsys):
        """Test that environment variable values with whitespace are handled correctly."""
        os.environ[self.test_envvar] = "  true  "
        
        result = util.get_envvar_switch(self.test_envvar, False, 'test setting')
        
        assert result is True
        captured = capsys.readouterr()
        assert f"ENVVAR {self.test_envvar} found, setting test setting to True." in captured.out
    
    def test_multiple_calls_independent(self):
        """Test that multiple calls to the function are independent."""
        os.environ[self.test_envvar] = 'true'
        result1 = util.get_envvar_switch(self.test_envvar, False)
        
        # Change environment variable
        os.environ[self.test_envvar] = 'false'
        result2 = util.get_envvar_switch(self.test_envvar, True)
        
        assert result1 is True
        assert result2 is False


    