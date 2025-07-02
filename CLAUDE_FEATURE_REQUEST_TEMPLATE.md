# Feature Request Template for Claude

## How to Use This Template

1. **Copy this file**: `cp CLAUDE_FEATURE_REQUEST_TEMPLATE.md my_feature_request.md`
2. **Fill in the sections** with your specific requirements
3. **Send the completed template** to Claude with: "Implement this feature request: [paste template content]"
4. **Or reference the file**: "Implement the feature described in my_feature_request.md"

Use this template to quickly request new features for the comparative genomics pipeline. Fill in the sections below to provide Claude with the context needed for accurate implementation.

## Basic Information

**Feature Name:** [Short descriptive name]

**Priority:** [High/Medium/Low]

**Component:** [Which part of the pipeline - e.g., visualization, client, service, analysis]

## Feature Description

**What:** [Brief description of what you want implemented]

**Why:** [Research goal or problem this solves]

**Where:** [Which files/modules should be modified or created]

## Technical Specifications

**Input Requirements:**
- Data format: [e.g., FASTA, JSON, CSV]
- Data source: [e.g., existing pipeline output, new API, user input]
- Dependencies: [any new libraries or external services needed]

**Output Requirements:**
- Format: [e.g., PNG plot, CSV file, JSON response]
- Location: [where output should be saved]
- Naming convention: [how files should be named]

**Processing Requirements:**
- Computational approach: [algorithm, statistical method, etc.]
- Performance considerations: [async/sync, memory usage, etc.]
- Error handling: [specific edge cases to consider]

## Integration Points

**Pipeline Step:** [Which of the 8 pipeline steps this fits into]
1. Ortholog Collection
2. Multiple Sequence Alignment  
3. Phylogenetic Tree Generation
4. Conservation Analysis
5. Variant Mapping
6. Scientific Visualization
7. Structure Retrieval
8. [New step]

**Data Dependencies:** [What existing pipeline outputs does this need]

**Configuration:** [Any new config parameters needed in genes_to_proteins.json or elsewhere]

## Scientific Context

**Research Question:** [How this relates to the core epilepsy gene research]

**Expected Results:** [What you expect to see/learn from this feature]

**Validation:** [How to verify the feature works correctly]

## Implementation Preferences

**Code Style:** [Any specific patterns to follow from existing codebase]

**Testing:** [Unit tests, integration tests, or manual testing needed]

**Documentation:** [README updates, plot annotations, etc.]

## Examples/References

**Similar Features:** [Point to existing code that does something similar]

**Data Samples:** [Provide example input/output if helpful]

**Research Papers:** [Any relevant citations or methods to follow]

## Quick Start Template

For simple requests, just fill this out:

```
Feature: [Name]
Component: [visualization/client/service/analysis]
Input: [what data]
Output: [what format, where]
Goal: [research objective]
Similar to: [existing feature in codebase]
```

## Usage Tips

1. **Be specific about file paths** - Use existing directory structure
2. **Reference existing patterns** - Point to similar code for consistency  
3. **Consider the research context** - Frame requests within epilepsy gene analysis
4. **Specify data flow** - How this connects to existing pipeline steps
5. **Include validation criteria** - How to know it worked correctly