# Python dockerfile image
FROM --platform=linux/amd64 python:3.9-slim

ENV PYTHONUNBUFFERED True

# Copy local code to the container image
ENV APP_HOME /app
WORKDIR $APPHOME
COPY . ./

# Install dependencies
RUN apt-get --allow-insecure-repositories update && \
    apt-get --assume-yes install gcc && \
    apt-get --assume-yes install pkg-config && \
    apt-get --assume-yes install libhdf5-serial-dev

RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Define gunicorn webserver exec parameters
CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 main:app