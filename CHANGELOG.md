# Changelog

## Develop
- Fix `mean_field_output_t` initialization, causing `start_time` to not be
  respected by the `user_stats` simulation component.
- Fix cyclic boundary rotation device bug, which tried to launch kernels
  with zero threads for ranks not containing cyclic boundaries.