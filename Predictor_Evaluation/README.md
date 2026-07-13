# This segment covers the relationship-prediction part of SNACKKSS. Specifically, re-running precision-recall metrics without smoothing ("unsmoothed"), testing our predictors against new manually curated relationship datasets ("new_gold_standard"), and testing DF1 against a Spearman correlation for signature-matching ("spearman")
# Go into the "import" directory, and gather the needed files specified in BUILD.sh

# Load the Docker image for SNACKKSS_Eval, then run:
docker run -v $(pwd)/:/app/ snackkss-eval-pipeline

