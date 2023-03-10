FROM tercen/runtime-flowsuite:3.15-2

COPY . /operator
WORKDIR /operator

ENV TERCEN_SERVICE_URI https://tercen.com

# RUN R -e "remotes::install_github('tercen/tim', force=TRUE)"
# RUN R -e "remotes::install_github('RGLab/FAUST', force=TRUE)"



ENTRYPOINT ["R", "--no-save", "--no-restore", "--no-environ", "--slave", "-f", "main.R", "--args"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]