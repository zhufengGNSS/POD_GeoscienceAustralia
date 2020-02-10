#!/usr/bin/env bash

# turn on strict bash
set -euo pipefail

#cd "$PROJECT_HOME"

organization=gamichaelmoore
repository=pod-pipelines

docker login -u "$DOCKER_HUB_USERNAME" -p "$DOCKER_HUB_PASSWORD"
#docker build -t "$organization/$repository" --build-arg AWS_SECRET_KEY="$AWS_SECRET_KEY_ID" --build-arg AWS_ACCESS_KEY="$AWS_ACCESS_KEY_ID" -f bitbucket-pipelines.dockerfile .
docker build -t "$organization/$repository" -f bitbucket-pipelines.dockerfile .
docker push "$organization/$repository"
