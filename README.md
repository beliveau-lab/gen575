# Ultrasound and Ultraviolet: Open Source Data

---

We have submitted our raw data (**UltrasoundandUltraviolet_Data.csv**), museum specimen information (**UltrasoundandUltraviolet_MuseumSpecimens.csv**), raw trimmed phylogenetic trees (**UltrasoundandUltraviolet_Trees.nex**), R script (**UltrasoundandUltraviolet_Script.R**), and photo deck of all ultraviolet-induced photoluminescence samples (including all species with no UVP; **UltrasoundandUltraviolet_UVP.pdf**).

## Descriptions

### **UltrasoundandUltraviolet**\_Data

*   *Name*: Currently valid scientific name
*   *Mass_g_Females*: Average mass of adult females (g)
*   *Mass_g_Males*: Average mass of adult males (g)
*   *Temporality*: The main temporal activity patten (Nocturnal or Diurnal)
*   *Gliding*: Whether or not this species can glide as a form of locomotion (Y/N)
*   *Openness*: Openness of the dominant habitat
    *   Closed = closed canopy, tall vegetation
    *   Open = open canopy, short or minimal vegetation
*   *Sociality*: Social interactions of adult individuals
    *   Social = Stable or cyclical social groups
    *   Solitary = Stable solitary living
*   *UV_Lit*: The dominant colour, if any, of ultraviolet-induced photoluminescence reported by other researchers in the literature; NA indicates that we could not find a record of ultraviolet-induced photoluminescence in the current literature.
*   *UV_Pub*: The dominant colour, if any, of ultraviolet-induced photoluminescence reported in this publication (photos available in **UltrasoundandUltraviolet_UVP**); NA indicates that we did not test for ultraviolet-induced photoluminescence in this publication.
*   *Lowest_kHz*: The absolute lowest frequency (kHz) of the dominant harmonic reported in literature (a list of all references and vocal ranges is provided in the supplemental information of the main article (Supplemental Information S1 and Supplemental References); NA indicates that we could not find vocal frequency range in the current literature.
*   *Highest_kHz*: The absolute highest frequency (kHz) of the dominant harmonic reported in literature (a list of all references and vocal ranges is provided in the supplemental information of the main article (Supplemental Information S1 and Supplemental References); NA indicates that we could not find vocal frequency range in the current literature.

### UltrasoundandUltraviolet\_MuseumSpecimens&#x20;

*   *Name*: Currently valid scientific name (note that this may not reflect the name used in the collections, we recommend using the ID to locate specimens in museum databases)
*   *Collection*: The museum that the specimen was stored in
*   *ID*: The catalogue ID for the specimen (this is the most reliable way to relocate specimens in both collections)
*   *Location of Origin*: The location that the specimen was obtained from (if available from the collection, NA indicates that location of origin was not catalogued by the collection). The city and province/state/territory were provided if possible, though some specimens were only listed with a country.
*   *Year of Origin*: The year that the specimen was obtained (if available from the collection, NA indicates that year of origin was not catalogued by the collection)
*   *Specimen Type*: The majority of specimens were pelts (which were sometimes stuffed), but we also note mounted specimens which were posed with additional features like glass eyes. Mounts are more likely to produce false positives as they are often treated with additional chemicals to achieve a desired look.
*   *Dorsal/Ventral_UV*: The conclusions reached by two independent observers as to the presence of ultraviolet-induced photoluminescence and the colour, if present. These conclusions are also present on each page of **UltrasoundUltraviolet_UVP**.

### **UltrasoundandUltraviolet**\_Trees

1000 node-dated trees were generated from the species list in **UltrasoundandUltraviolet_Data**. Some species had nomenclature inconsistent with our dataset:

*   *Lagothrix lagothricha*, Vertlife = *Lagothrix lagotricha*
*   *Micoureous demerarae*, Vertlife = *Marmosa demerarae*
*   *Paragalago* spp., Vertlife = *Galagoides* spp.
*   *Plecturocebus moloch*, Vertlife = *Callicebus moloch*

Additionally, *Petaurus notatus* is a recently described species that was not available at Vertlife at the time of publication; this species was added at 1MYA to the *Petaurus breviceps* tip using **UltrasoundUltraviolet_Script**.

### **UltrasoundandUltraviolet**\_UVP

Museum specimens photographed under normal conditions (white light), ultraviolet light, ultraviolet light + yellow gel filter (to reduce the amount of excess purple-blue light), and the filtered photograph with a slight edit to reduce the excess of yellow (details provided in the supplemental information of the main article (SupplementalInformationS1)). The ventral and dorsal sides were photographed for all specimens; museum collection indicated in the bottom-right corner.

## Key Information Sources

Body Mass, Habitat Openness, Sociality, Temporality, and Gliding were derived from the following sources:

*   Animal Diversity Web
*   Integrative Taxonomic Information System
*   IUCN Red List

Phylogenetic trees were trimmed from the Mammalian Super Tree:

*   Vertlife

## Code/Software

R is required to run *Nocturnal Gliders_Script*; the script was created using version 4.2.3.
Annotations are provided throughout the script through 1) library loading, 2) dataset loading and cleaning, 3) analyses, and 4) figure creation.
