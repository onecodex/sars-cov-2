test:
	py.test -vv tests/
	@echo "Successfully ran all tests."

lint:
	pre-commit run --all-files
	@echo "Successfully linted all files."
