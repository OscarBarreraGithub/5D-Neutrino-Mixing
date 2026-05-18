import { defineCollection, z } from 'astro:content';
import { glob } from 'astro/loaders';

// ---------------------------------------------------------------------------
// Schema for one ingested catalog entry.  Mirrors the YAML sidecar fields the
// ingestion script (scripts/ingest_catalog.py) hoists out of
// flavor_catalog/processes/**/*.yaml plus the .tex prose sections.  Every
// field is optional except `process_id` and `family` because the catalog is
// heterogeneous across waves.
// ---------------------------------------------------------------------------

const pdgValue = z
  .object({
    label: z.string().nullish(),
    observable: z.string().nullish(),
    value: z.union([z.string(), z.number()]).nullish(),
    uncertainty: z.union([z.string(), z.number()]).nullish(),
    units: z.string().nullish(),
    display: z.string().nullish(),
    year: z.union([z.string(), z.number()]).nullish(),
    experiment: z.string().nullish(),
    source_url: z.string().nullish(),
    snapshot_path: z.string().nullish(),
    sha256: z.string().nullish(),
    access_date: z.string().nullish(),
    cl: z.string().nullish(),
    limit_type: z.string().nullish(),
  })
  .passthrough();

const sectionsSchema = z
  .object({
    process: z.string().nullish(),
    pdg_prose: z.string().nullish(),
    relevance: z.string().nullish(),
    post_2008: z.string().nullish(),
    validity: z.string().nullish(),
    code_coverage_prose: z.string().nullish(),
    implementation_difficulty_prose: z.string().nullish(),
    key_references: z.string().nullish(),
  })
  .passthrough();

const codeCoverageSchema = z
  .object({
    status: z.enum(['YES', 'PARTIAL', 'NO']).or(z.string()).nullish(),
    coverage_alias: z.string().nullish(),
    searched_directories: z.array(z.unknown()).nullish(),
    grep_commands: z.array(z.unknown()).nullish(),
    evidence: z.array(z.unknown()).nullish(),
  })
  .passthrough();

const entrySchema = z
  .object({
    process_id: z.string(),
    family: z.string(),
    family_original: z.string().nullish(),
    tier: z.enum(['PRIMARY', 'SECONDARY']).default('PRIMARY'),
    standard_notation: z.string().nullish(),
    process_name: z.string().nullish(),
    owner: z.string().nullish(),
    fact_check_verdict: z
      .enum(['VERIFIED', 'PARTIAL', 'MISMATCH', 'FAILED', 'UNKNOWN'])
      .or(z.string())
      .nullish(),
    implementation_difficulty: z.string().nullish(),
    implementation_difficulty_reason: z.string().nullish(),
    cycle_count: z.number().int().nullish(),
    priority_rationale: z.string().nullish(),
    promoted_in_wave: z.union([z.string(), z.number()]).nullish(),
    code_coverage: codeCoverageSchema.nullish(),
    code_coverage_status: z.string().nullish(),
    pdg_values: z.array(pdgValue).default([]),
    supporting_measurements: z.unknown().nullish(),
    paper_era_reference: z.unknown().nullish(),
    theory_context: z.unknown().nullish(),
    auxiliary_theory_inputs: z.unknown().nullish(),
    post_2010_developments: z.unknown().nullish(),
    source_shas: z.record(z.string()).nullish(),
    source_shas_count: z.number().int().default(0),
    access_dates: z.array(z.string()).default([]),
    latest_status: z.string().nullish(),
    status_history: z.unknown().nullish(),
    sections: sectionsSchema.default({}),
    source_yaml: z.string().nullish(),
    source_tex: z.string().nullish(),
    worklog_path: z.string().nullish(),
  })
  .passthrough();

const entries = defineCollection({
  loader: glob({ pattern: '*.json', base: './src/content/entries' }),
  schema: entrySchema,
});

export const collections = { entries };
