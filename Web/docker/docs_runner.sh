set -v
set -e

cat build_docs.sh

# Copy the REFPROP files here
cp -r ${HOME}/REFPROP_sources .

# Run the build of the docs
docker-compose build
docker-compose run worker bash /coolprop/Web/docker/build_docs.sh $1
