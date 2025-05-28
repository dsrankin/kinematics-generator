## Kinematic Generator

Running the test:
<pre>
root -l test.cc
</pre>

#### Specifying decays

The current script will always start with the production of and A and B.
You can then specify the decay chain through the maps at the top of `test.cc`

This specifies if the particle is stable or not.
```
std::map<char, bool> pstable = {
    { 'A', false },
...
```

This specifies if the particle is visible or not. For unstable particles this has no effect, but for stable particles it will result in the particle being included in the MET.
```
std::map<char, bool> pvisible = {
    { 'A', false },
...
```

This specifies what the particle should decay into.
```
std::map<char, std::pair<char,char>> pdecay = {
    { 'A', {'W', 'b'} },
...
```

This specifies the mass of the particle. 
```
std::map<char, double> pmass = {
    { 'A', 5.5 },
...
```

Note: a three-body decay can be done by specifying a decay such that the daughter masses are larger than that of the mother.
In this case, the heavier mass daughter will immediate be decayed.
For example, if A (m=100) decays into C (m=50) and D (m=500), and D decays into E (m=5) and F (m=0), then this will actuall perform a 3 body decay of A into C, E, and F.

This specifies the ID number to use to encode the particle.
```
std::map<char, int> pid = {
    { 'A', 9999 },
...
```
Finally, the number of events can be specifies with this line: `int nevents = 1000000;`

For more, better just read the code.
