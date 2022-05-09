
# build

```shell
VERSION=0.14.0
docker build -t tercen/simple_docker_operator:$VERSION .
docker push tercen/simple_docker_operator:$VERSION
# update operator.json file with correct docker image version
{
        echo '{'
        echo '"name": "simple docker operator",'
        echo '"description": "simple docker operator",'
        echo '"tags": [""],'
        echo '"authors": ["tercen"],'
        echo '"urls": ["https://github.com/tercen/simple_docker_operator"],'
        echo '"container":"tercen/simple_docker_operator:'$VERSION'",'  
        echo '"properties": [ ]'
        echo '}'
} > operator.json

git add -A && git commit -m "$VERSION" && git tag $VERSION && git push && git push --tags

```

# inspect

```shell
docker run --rm --entrypoint=bash tercen/simple_docker_operator:$VERSION -c "R --version"
docker run -it --rm --entrypoint=bash tercen/simple_docker_operator:$VERSION
```
 
# push

```shell
docker push tercen/simple_docker_operator:$VERSION
```