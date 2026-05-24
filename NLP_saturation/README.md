#Go into the "import" directory, and gather the needed files specified in BUILD.sh

#raw/ML_steps.txt is a reference file, specifying the number of steps taken by each training run in this study. We calculate the number of steps per checkpoint by dividing this number by six and rounding down to the integer, for each task/database pair (e.g. CREEDS, study classification), thus ensuring that there will always be at least five checkpoints.

#Load the Docker image for SNACKKSS_NLP, then run:
docker run -v $(pwd)/:/app/ snackkss-nlp-pipeline

