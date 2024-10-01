# Integration tests

Small cases, which tests that overall functionality have not been broken.
Ideally each test should not be longer than a few seconds.
The tests are run with pytest, you can install it using pip.

Execute
```
pytest --tb=no
```
to run the tests.


Directory structure.

- meshes. Mesh files that can be used in the tests. Add a new one if needed.
- case\_tempates. A few case files that can be uses as a start to build one for your case.
- Each test suite should be in its own directory, where you can also add necessary stuff such as reference log files.
- Each test suite should be in a file test_\*.py, this way it is automatically discovered by pytest.
- Each test in the suite should be in a function which starts with `test_`.
