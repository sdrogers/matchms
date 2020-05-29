```
$ git clone https://github.com/matchms/matchms.git
$ git checkout <some sha>
$ SHA=`git log --format=%h -n 1`
$ cd matchms
$ conda env create --file conda/environment-dev.yml
$ conda activate matchms-dev
$ pip install pyprof2calltree
$ pip install --editable .
$ cd profiling
$ python -m cProfile -o pprof.user_workflow.$SHA profile_user_workflow.py
$ pyprof2calltree -i pprof.user_workflow.$SHA -o cachegrind.out.user_workflow.$SHA
$ sudo apt install kcachegrind

```

Use ``kcachegrind`` to open ``cachegrind.out.user_workflow.$SHA``.
