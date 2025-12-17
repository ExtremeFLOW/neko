# Changelog

## Develop
- Error and warning routines are now hooked to pFUnit's exceptions, making it
  possible to test for errror emission.
- Fix `mean_field_output_t` initialization, causing `start_time` to not be
  respected by the `user_stats` simulation component.
