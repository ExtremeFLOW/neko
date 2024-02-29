# Immersed boundary of a Stanford Bunny

This example illustrates the brinkman condition source term used to simulate an
immersed boundary. A Stanford Bunny model is placed in side a square duct and
flow is simulated.

```
Bunny url:              https://www.thingiverse.com/thing:151081
Bunny domain:           [ [-24, 85], [-41, 45], [5, 112] ]
Design domain:          [ [-200, 600], [-100, 100], [0, 200] ]
```

The run script will construct a design domain as above with 16 x 8 x 8 elements.
The user can decide to specify another resolution from the terminal.

The user can change how sharp the interface between fluid and solid is by
changing the `distance_transform.value`, the value specifies the width of the
transition region from no to full resistance. The size of the Brinkman
resistance limits should be adjusted if the Reynolds number is modified.
