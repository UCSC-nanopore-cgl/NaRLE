# Definitions
build_output = runtime/rtg-tools.jar
runtime_fullpath = $(realpath runtime)
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = tpesout/rtg_tools
tag = 3.8.3--${git_commit}

build: .PHONY
	# build with all that java trash
	cd build && docker build -t rtgtoolsbuild .
	docker run -v ${runtime_fullpath}:/data rtgtoolsbuild cp /home/rtgtools-protected/rtg-tools.jar /home/rtgtools-protected/LICENSE.txt /data
	# build with just the jar
	cd runtime && docker build -t ${name}:${tag} .
	-docker rmi ${name}:latest
	docker tag ${name}:${tag} ${name}:latest
	-docker rmi -f rtgtoolsbuild

push: build
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest

.PHONY:
