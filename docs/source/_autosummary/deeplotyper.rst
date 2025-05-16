deeplotyper
===========

.. automodule:: deeplotyper

   
   .. rubric:: Functions

   .. autosummary::
   
      TypedDict
      apply_alignment_gaps
      build_linear_coords
      build_raw_genome_coords
      build_raw_transcript_coords
      dataclass
      field
      find_orfs
      get_longest_orf
      make_aligner
   
   .. rubric:: Classes

   .. autosummary::
   
      Base
      BaseCoordinateMapping
      CodonCoordinateMapping
      CodonIndex
      ExonDef
      HaplotypeEvent
      HaplotypeRemapper
      MappingEntry
      NewTranscriptSequences
      Position
      RawBase
      Region
      SequenceCoordinateMapper
      SequenceMappingResult
      TranscriptMappingResult
   
.. rubric:: Modules

.. autosummary::
   :toctree:
   :recursive:

   alignment_utils
   coordinate_utils
   data_models
   extras
   mapper
   orf_utils
   remapper
   vcf_haploevents
