## Kinematic Generator

Running the test:
<pre>
root -l test.cc
</pre>
or with an explicit process card:
<pre>
root -l -q 'test.cc("process.card")'
</pre>

#### Process card configuration

The generator now reads all model/process settings from a process card
(`process.card` by default), so you do not need to edit `test.cc`.

Card syntax:
```
nevents <int>
outfile <root_file_name>
production <particleA> <particleB>
particle <name> <mass> <pid> <stable(0|1)> <visible(0|1)> <decay1> <decay2>
```

For stable particles, `decay1` and `decay2` are ignored, but placeholders are required.
See `process.card` for a full example.
