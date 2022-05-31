FROM tercen/runtime-flowsuite:3.15-1

RUN R -e "remotes::install_github('RGLab/scamp', ref='v0.5.0')"
RUN R -e "remotes::install_github('RGLab/FAUST', ref='f900614')"

COPY . /operator
WORKDIR /operator

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT ["R", "--no-save", "--no-restore", "--no-environ", "--slave", "-f", "main.R", "--args"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]