# Definitions
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
version = 2.0.5
name = tpesout/alignqc
tag = ${version}--${git_commit}

build:
	# do the docker build
	docker build --build-arg git_commit=${git_commit} -t ${name}:${tag} .
	docker tag ${name}:${tag} ${name}:latest

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest
