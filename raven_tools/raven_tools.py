import raven_model as rm

suffix = [
    "rvi", "rvh", "rvp", "rvc", "rvt"
]

models = [
    "GR4J",
    "HYMOD",
    "HMETS",
    "HBV",
    "MOHYSE"
]
for m in models:
    model_instance = rm.RavenModel(model_type=m, catchment="Broye")
    model_instance.create_dirs()
    for s in suffix:
        model_instance.write_rvx(ostrich_template=True, rvx_type=s)
    model_instance.write_ost()

import os
src = '/home/mainman/PycharmProjects/raven-tools/raven_tools/preappend_dblquotes.py'
dst = '/home/mainman/PycharmProjects/raven-tools/raven_tools/preappend_dblquotes_symlink.py'
os.symlink(src, dst)