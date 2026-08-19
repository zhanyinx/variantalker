"""
Test runner for MAFigate test suite.
Run all tests or specific test modules.
"""

import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_all_tests():
    """Run all tests in the test suite."""
    print("🧪 Running MAFigate Test Suite")
    print("=" * 50)

    # Discover and run all tests
    loader = unittest.TestLoader()
    start_dir = os.path.dirname(os.path.abspath(__file__))
    suite = loader.discover(start_dir, pattern="test_*.py")

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print("\n" + "=" * 50)
    print(f"📊 Test Summary:")
    print(f"   Tests run: {result.testsRun}")
    print(f"   Failures: {len(result.failures)}")
    print(f"   Errors: {len(result.errors)}")
    print(f"   Skipped: {len(result.skipped)}")

    if result.failures:
        print(f"\n❌ Failures:")
        for test, traceback in result.failures:
            print(f"   - {test}: {traceback.split('AssertionError:')[-1].strip()}")

    if result.errors:
        print(f"\n💥 Errors:")
        for test, traceback in result.errors:
            print(f"   - {test}: {traceback.split('Exception:')[-1].strip()}")

    success = len(result.failures) == 0 and len(result.errors) == 0
    print(f"\n{'✅ All tests passed!' if success else '❌ Some tests failed!'}")

    return success


def run_specific_test(test_module):
    """Run a specific test module."""
    print(f"🧪 Running {test_module}")
    print("=" * 50)

    try:
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromName(test_module)
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)

        success = len(result.failures) == 0 and len(result.errors) == 0
        print(f"\n{'✅ Tests passed!' if success else '❌ Tests failed!'}")
        return success

    except Exception as e:
        print(f"❌ Error running {test_module}: {e}")
        return False


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run specific test module
        test_module = sys.argv[1]
        if not test_module.startswith("test_"):
            test_module = f"test_{test_module}"
        run_specific_test(test_module)
    else:
        # Run all tests
        run_all_tests()
