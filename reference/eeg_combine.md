# Combine `eegUtils` objects

Combine multiple `eeg_epochs`, `eeg_data`, or `eeg_evoked` objects into
a single object. The function will check the `participant_id` entry in
the `epochs` structure of each object to see if the objects come from a
single participant or from multiple participants. If the data are from
multiple participants, it will create an `eeg_group` object. For
individual participants, it will check for duplicate epochs. For most
objects, it will concatenate the objects if none are found. However, for
`eeg_data` it will instead try to correct the epoch numbers. Check the
details below for further advice.

## Usage

``` r
eeg_combine(data, ...)

# S3 method for class 'list'
eeg_combine(data, ...)

# S3 method for class 'eeg_data'
eeg_combine(data, ..., check_timings = TRUE)

# S3 method for class 'eeg_epochs'
eeg_combine(data, ..., check_timings = TRUE)

# S3 method for class 'eeg_evoked'
eeg_combine(data, ...)
```

## Arguments

- data:

  An `eeg_data`, `eeg_epochs`, or `eeg_evoked` object, or a list of such
  objects.

- ...:

  additional `eeg_data` or `eeg_epochs` objects

- check_timings:

  Check whether sample times / epoch numbers are continuously ascending;
  if not, modify so that they are. Useful when, for example, combining
  epochs derived from multiple recording blocks. Defaults to TRUE

## Value

If all objects have the same `participant_id`, returns an object of the
same class as the original input object. If the objects have different
`participant_id` numbers, an object of both class `eeg_group` and the
same class as the original input object.

## Methods (by class)

- `eeg_combine(list)`: Method for combining lists of `eeg_data` and
  `eeg_epochs` objects.

- `eeg_combine(eeg_data)`: Method for combining `eeg_data` objects.

- `eeg_combine(eeg_epochs)`: Method for combining `eeg_epochs` objects

- `eeg_combine(eeg_evoked)`: Method for combining `eeg_evoked` objects

## Combining `eeg_data` objects

Combining `eeg_data` is mainly intended to be used for combining
multiple recordings from a single participant prior to subsequent
epoching. Thus, `check_timings` defaults to true, and the function will
change the epochs and timing structures of the resulting combined object
to be as if it were a single recording. The objects will be combined in
the input order, so ensure that the objects are input in chronological
order.

## Combining `eeg_epochs` objects

There are several scenarios where you might wish to combine
`eeg_epochs`. For example, a user may have processed continuous data in
smaller chunks reflecting short recording blocks before epoching. They
then wish to combine these into a single object. In that case, the epoch
numbering should reflect chronological ordering and needs to be
corrected.

If `check_timings == TRUE`, the function will perform several checks
before combining objects. First, it will check for duplicate epochs in
the `epochs` structure of each object. If each object only has unique
epochs, the objects will be combined without correction. Thus, combining
across separate recordings or separate participants will not elicit
correction. The user should ensure

If there are any duplicates (e.g. a participant has more than one epoch
numbered one from the same recording), it will then check if there are
any missing epochs. If there are, the new trial numbering cannot be
automatically determined, so the objects cannot be combined without
further manual intervention. If there are no missing epochs, it will
then check if there is any decreases in epoch numbers across objects. If
there are any, then the epoch numbers and timings for the objects will
be adjusted.

Alternatively, the user may wish to combine `eeg_epochs` objects from
different participants or from entirely different recording sessions of
the same participant. In this case, no correction of timings or epoch
numbers is desirable. `check_timings == TRUE` should detect this and
skip correction, but can be explicitly set to `FALSE`.

## Author

Matt Craddock, <matt@mattcraddock.com>
