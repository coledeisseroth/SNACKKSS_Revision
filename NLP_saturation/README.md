#Go into the "import" directory, and gather the needed files specified in BUILD.sh

#raw/ML_steps.txt is a reference file, specifying the number of steps taken by each training run in this study.

#Load the Docker image for SNACKKSS_NLP, then run:
docker run -v $(pwd)/:/app/ snackkss-nlp-pipeline

