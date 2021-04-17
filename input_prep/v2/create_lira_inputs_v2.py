#!/usr/bin/env python

from TaskRunner import TaskRunner

runner = TaskRunner('lira_input_prep.yaml')
runner.create_LIRA_inputs()