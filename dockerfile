FROM continuumio/anaconda3
COPY ./seir_interactive.py /root/
CMD bokeh serve /root/seir_interactive.py --port 80