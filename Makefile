test:
	py.test -vv tests/
	@echo "Successfully ran all tests."

lint:
	pre-commit run --all-files
	@echo "Successfully linted all files."

build: build/local build/docker
	@echo "Set up environment and built Docker image"

build/local:
	pip install -r requirements.txt

build/docker:
	docker build -t covid19 .
