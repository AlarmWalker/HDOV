# Spatial-Temporal Index for Immersive Virtual 3D Spaces

## Overview

This project introduces a spatial-temporal index tailored for immersive virtual 3D environments, such as those found in the Metaverse. The Java-based software application manages the spatial index of a three-dimensional urban model, focusing on Den Haag. It aims to enhance user interactions and scalability within these environments.

## Authors

- Patrick Kim ([gkim2@unb.ca](mailto:gkim2@unb.ca))
- Sadhana Suresh Chettiar ([sadhana.chettiar@unb.ca](mailto:sadhana.chettiar@unb.ca))
- Shidur Sharma Durba ([shidursharma.durba@unb.ca](mailto:shidursharma.durba@unb.ca))

## Abstract

The advent of the Metaverse has sparked a pressing need for efficient management of entities within immersive virtual 3D environments. This study introduces a pioneering approach to implementing a spatial-temporal index tailored for such environments, aiming to enhance user interactions and scalability. Leveraging the unique characteristics of immersive virtual spaces, our method showcases a Java-based software application designed to manage the spatial index of a three-dimensional urban model, focusing on Den Haag.

## Key Features

- **Spatial Indexing:** Utilizes STRtree and BBx-index structures for efficient spatial indexing.
- **Temporal Indexing:** Implements temporal indexing to support queries regarding past, present, and future positions of objects.
- **Database Integration:** Uses PostgreSQL with the PostGIS extension and 3DCityDB for managing 3D urban models.
- **Visibility Analysis:** Calculates the Degree of Visibility (DoV) of objects within the urban model from various viewpoints.

## Technologies Used

- Java
- PostgreSQL with PostGIS
- 3DCityDB
- JTS (Java Topology Suite)
- JSON/CityJSON

## Project Setup

1. **Database Creation:**
   - Establish a PostgreSQL database with PostGIS extension.
   - Populate the schema using the CREATE_DB.sh script provided by 3DCityDB.
   - Import the LOD2 dataset using the Importer/Exporter tool.

2. **Data Retrieval:**
   - Extract geometry objects and associated data from the PostgreSQL database.
   - Convert POLYHEDRALSURFACE Z geometries to LINESTRING Z and POLYGON Z using the ST_ConvexHull function.

3. **Implementation:**
   - Develop a Java application to manage the spatial index.
   - Implement the STRtree and BBx-index structures for spatial indexing.
   - Conduct visibility analysis using grid cells and line-of-sight vectors.

## Visibility Calculation Methods

- **STRtree Method:** Utilizes a spatial index to efficiently query nearby geometries and perform visibility calculations.
- **Brute Force Method:** Iterates through all geometries without using a spatial index, providing a baseline for comparison.

## Evaluation

- Conducted experiments using the DenHaag dataset containing over 100,000 building models.
- Compared the performance of the STRtree-based approach with a brute-force approach.
- Measured execution time, number of visible buildings, and overall efficiency.

## Results

- The STRtree-based method significantly outperformed the brute-force approach in terms of processing time and efficiency.
- Implemented a temporal index using the BBx-tree structure to handle dynamic changes in spatial data over time.

## Future Work

- Optimize the temporal indexing algorithm.
- Explore alternative techniques for managing spatial indices within 3D immersive environments.

## References

- [3D City Database Documentation](https://3dcitydb-docs.readthedocs.io/en/latest/index.html)
- Lin, D., Jensen, C., Ooi, B., & Šaltenis, S. (2005). Efficient indexing of the historical, present, and future positions of moving objects. *Proceedings of the ACM International Conference on Management of Data (SIGMOD)*.
- Shou, L., Huang, Z., & Tan, K. L. (2003). HDoV-tree: the structure, the storage, the speed. *Proceedings of the IEEE International Conference on Data Engineering (ICDE)*.

---

**For more details, please refer to the [project report](./spatial-temporal index.pdf).**
