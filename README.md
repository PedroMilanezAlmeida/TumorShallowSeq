Cancer prognosis with shallow tumor RNA sequencing
==================================================

Here you will find out how to perform the analysis at the core of our
2020 Nature Medicine paper:
<a href="https://www.nature.com/articles/s41591-019-0729-3" class="uri">https://www.nature.com/articles/s41591-019-0729-3</a>.

Your R version has to be 3.6.1 or higher for this to work.

    if(sessionInfo()$R.version$version.string < '3.6.1'){
      stop(paste0("\n##### This will not run for R versions older than 3.6.1 #####\n",
                  "##### your version is: ",
                  sessionInfo()$R.version$version.string,
                  " #####\n",
                  "##### update R and rerun"))
    }
