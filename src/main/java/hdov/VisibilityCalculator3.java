package hdov;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import java.sql.*;
import java.time.LocalDate;
import java.time.Month;
import java.util.*;
import java.util.Date;

class GeometryWithDate {
    private Geometry geometry;
    private Date date;

    public GeometryWithDate(Geometry geometry, Date date) {
        this.geometry = geometry;
        this.date = date;
    }

    public Geometry getGeometry() {
        return geometry;
    }

    public Date getDate() {
        return date;
    }
}

// Vertical Storage Scheme code starts
class VPage {
    // Represents a V-page containing visibility data
    private Map<Integer, Double> visibilityData; // Maps object IDs to their DoV values

    public VPage() {
        this.visibilityData = new HashMap<>();
    }

    public void addVisibilityData(int objectId, double dov) {
        visibilityData.put(objectId, dov);
    }

    public Double getVisibility(int objectId) {
        return visibilityData.getOrDefault(objectId, 0.0);
    }

    public Map<Integer, Double> getVisibilityData() {
        return visibilityData;
    }
}

class VPageIndex {
    // Represents the V-page-index
    private List<VPage> vPages;

    public VPageIndex(int size) {
        vPages = new ArrayList<>(Collections.nCopies(size, null));
    }

    public void setVPage(int index, VPage vPage) {
        vPages.set(index, vPage);
    }

    public VPage getVPage(int index) {
        return vPages.get(index);
    }
    public int getSize() {
        return vPages.size();
    }
    public List<VPage> getVPages() {
        return vPages;
    }

    //Addition of code to address problem of "Reached the limit of VPageIndex. Some nodes will not be processed" starts
    public void resize(int newSize) {
        if (newSize > vPages.size()) {
            vPages.addAll(Collections.nCopies(newSize - vPages.size(), null));
        }
    }
    //Addition of code to address problem of "Reached the limit of VPageIndex. Some nodes will not be processed" ends

}

//Vertical Storage scheme code ends

//HDoV Logical Structure code starts
class HDoVTreeNode {
    double aggregatedDoV; // Aggregated DoV for internal nodes, individual DoV for leaf nodes
    Geometry mbr; // Minimum bounding rectangle
    List<HDoVTreeNode> children; // Child nodes
    double lod; // Level of Detail

    //HDoV search algorithm code starts
    int polygonCount; // Field to store the polygon count
    //HDoV search algorithm code ends

    Date creationDate;

    Date deletionDate;

    // Vertical Storage Scheme code starts
    private int vPageIndexOffset; // Offset in the V-page-index
    // Method to set the V-page-index offset

    private int id; // ID field

    public void setVPageIndexOffset(int offset) {
        this.vPageIndexOffset = offset;
    }

    // Method to get the V-page-index offset
    public int getVPageIndexOffset() {
        return vPageIndexOffset;
    }

    // Getter method for the ID
    public int getId() {
        return id;
    }

    // Setter method for the ID (if needed)
    public void setId(int id) {
        this.id = id;
    }

    // Vertical Storage Scheme code ends

    // Constructor for leaf nodes
    public HDoVTreeNode(double dov, Geometry geometry, double lod, int polygonCount, Date creationDate) //note: int polygonCount added for search algorithm rest others were already part of HDoV logical structure
    {
        this.aggregatedDoV = dov;
        this.mbr = geometry;
        this.lod = lod;
        this.children = new ArrayList<>();
        //HDoV search algorithm code starts
        this.polygonCount = polygonCount;
        //HDoV search algorithm code ends
        this.creationDate = creationDate;
        this.deletionDate = new Date();
    }

    // Constructor for internal nodes
    public HDoVTreeNode(List<HDoVTreeNode> children) {
        this.children = children;
        this.aggregatedDoV = calculateAggregatedDoV(children);
        // Calculate MBR for internal node...
        this.mbr = calculateMBR(children);
        this.lod = 0; // LoD is not relevant for internal nodes
    }

    private double calculateAggregatedDoV(List<HDoVTreeNode> children) {
        double sum = 0;
        for (HDoVTreeNode child : children) {
            sum += child.aggregatedDoV;
        }
        return sum;
    }

    private Geometry calculateMBR(List<HDoVTreeNode> children) {
        GeometryFactory geometryFactory = new GeometryFactory();
        if (children.isEmpty()) {
            //return geometryFactory.createPolygon(); // Return an empty polygon if there are no children // Original code
            //Test code
            return VisibilityCalculator3.geometryFactory.createPolygon();
        }

        // Initialize min and max coordinates
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        double maxY = Double.MIN_VALUE;

        for (HDoVTreeNode child : children) {
            Envelope env = child.mbr.getEnvelopeInternal();
            minX = Math.min(minX, env.getMinX());
            minY = Math.min(minY, env.getMinY());
            maxX = Math.max(maxX, env.getMaxX());
            maxY = Math.max(maxY, env.getMaxY());
        }

        // Create a polygon that represents the MBR encompassing all child MBRs
        Coordinate[] coordinates = new Coordinate[]{
                new Coordinate(minX, minY),
                new Coordinate(maxX, minY),
                new Coordinate(maxX, maxY),
                new Coordinate(minX, maxY),
                new Coordinate(minX, minY) // Closed linear ring
        };

        return geometryFactory.createPolygon(coordinates);
    }

    // Method to cluster leaf nodes based on proximity
    private static List<List<HDoVTreeNode>> clusterLeafNodes(List<HDoVTreeNode> leafNodes) {
        // Placeholder for clustering logic: this should be replaced with actual clustering logic
        // For now, we split the list into smaller lists (clusters) of equal size for simplicity
        List<List<HDoVTreeNode>> clusters = new ArrayList<>();
        int clusterSize = 10; // Adjust cluster size as needed
        for (int i = 0; i < leafNodes.size(); i += clusterSize) {
            clusters.add(leafNodes.subList(i, Math.min(i + clusterSize, leafNodes.size())));
        }
        return clusters;
    }

    // Method to create internal nodes from clusters and build the HDoV-tree
    public static HDoVTreeNode buildHDoVTree(List<HDoVTreeNode> leafNodes) {
        List<List<HDoVTreeNode>> clusters = clusterLeafNodes(leafNodes);
        List<HDoVTreeNode> internalNodes = new ArrayList<>();

        for (List<HDoVTreeNode> cluster : clusters) {
            internalNodes.add(new HDoVTreeNode(cluster));
        }

        // Recursively cluster internal nodes until one root node remains
        while (internalNodes.size() > 1) {
            clusters = clusterLeafNodes(internalNodes);
            internalNodes = new ArrayList<>();
            for (List<HDoVTreeNode> cluster : clusters) {
                internalNodes.add(new HDoVTreeNode(cluster));
            }
        }

        return internalNodes.isEmpty() ? null : internalNodes.get(0);
    }

    public void traverseAndPrint(Coordinate viewpoint) {
        // Priority based on proximity to viewpoint
        double distance = this.mbr.distance(new GeometryFactory().createPoint(viewpoint));
        System.out.println("Node distance from viewpoint: " + distance + ", DoV: " + aggregatedDoV);

        // Recursively traverse children sorted by proximity
        children.stream()
                .sorted(Comparator.comparingDouble(child -> child.mbr.distance(new GeometryFactory().createPoint(viewpoint))))
                .forEach(child -> child.traverseAndPrint(viewpoint));
    }

    public void printNode() {
        System.out.println("Node DoV: " + aggregatedDoV + ", MBR: " + mbr.toText() + ", LoD: " + lod + ", Date: " + creationDate);
        for (HDoVTreeNode child : children) {
            child.printNode(); // Recursive call for child nodes
        }
    }

    //HDoV search algorithm code starts
    // Method to search the HDoV tree
    public void searchHDoVTree(Coordinate queryPoint, List<GeometryWithDate> results, double eta) {
        if (this.aggregatedDoV == 0) return; // Check if DoV is 0

        if (this.children.isEmpty()) { // Leaf node
            results.add(new GeometryWithDate(this.mbr, this.creationDate)); // Add MBR of the leaf node to results
        } else { // Internal node
            double nodeDoV = this.aggregatedDoV;
            int nvo = this.children.size(); // Number of visible objects (children)

            // Calculate 's' and 'm' based on your data structure and requirements
            double s = calculatePolygonRatio(); // Implement this method based on your data
            int m = countLeafDescendants(); // Implement this method to count leaf descendants

            double h = Math.log(m) / Math.log(2); // Height of the subtree

            if (nodeDoV <= eta && (h * (1 + Math.log(s)) < Math.log(nvo))) {
                results.add(new GeometryWithDate(this.mbr, this.creationDate)); // Add MBR of the internal node to results
            } else {
//                for (HDoVTreeNode child : this.children) {
//                    double dov = calculateDoV(child, queryPoint); // Implement this method to calculate DoV based on queryPoint
//
//                    // Check if the node's DoV is greater than the threshold eta
//                    if (dov > eta) {
//                        // Logic for nodes with high DoV - retrieve or traverse with high detail
//                        child.searchHDoVTree(queryPoint, results, eta);
//                    } else {
//                        // Logic for nodes with low DoV - retrieve low-level internal LoD
//                        // Add child's low-detail representation to results
//                    }
//                }

                for (HDoVTreeNode child : this.children) {
                    child.searchHDoVTree(queryPoint, results, eta); // Recurse into children
                }
            }
        }
    }
    //HDoV search algorithm code ends

    // Vertical Storage Scheme code starts

//    public void searchHDoVTree(Coordinate queryPoint, List<Geometry> results, double eta, VPageIndex vPageIndex) {
//        VPage vPage = vPageIndex.getVPage(this.vPageIndexOffset);
//        if (vPage == null || this.aggregatedDoV == 0) return; // Check visibility from V-page
//
//        // ... [rest of the search logic] ...
//    }

    public void searchHDoVTree(Coordinate queryPoint, List<GeometryWithDate> results, double eta, VPageIndex vPageIndex) {
        VPage vPage = vPageIndex.getVPage(this.vPageIndexOffset);
        if (vPage == null) return; // If no V-page, the branch is not visible

        // Retrieve the DoV value for this node from the V-page
        double nodeDoV = vPage.getVisibility(this.getId()); // Assuming each node has a unique ID

        if (this.children.isEmpty()) { // Leaf node
            if (nodeDoV > 0) {
                // Add MBR of the leaf node to results
                results.add(new GeometryWithDate(this.mbr, this.creationDate));
            }
        } else { // Internal node
            int nvo = countVisibleObjects(nodeDoV, this.children); // Number of visible objects (children)

            double s = calculatePolygonRatio(); // Polygon ratio
            int m = countLeafDescendants(); // Count of leaf descendants
            double h = Math.log(m) / Math.log(2); // Height of the subtree

            if (nodeDoV <= eta && (h * (1 + Math.log(s)) < Math.log(nvo))) {
                results.add(new GeometryWithDate(this.mbr, this.creationDate)); // Add MBR of the internal node to results
            } else {
                // Recurse into children
                for (HDoVTreeNode child : this.children) {
                    child.searchHDoVTree(queryPoint, results, eta, vPageIndex);
                }
            }
        }
    }

    // Helper method to count the number of visible objects based on the DoV value
    private int countVisibleObjects(double nodeDoV, List<HDoVTreeNode> children) {
        int count = 0;
        for (HDoVTreeNode child : children) {
            if (child.aggregatedDoV > nodeDoV) {
                count++;
            }
        }
        return count;
    }

    // Vertical Storage Scheme code ends

    //HDoV search algorithm code starts
    // Calculate the ratio of polygons in this node to the sum of polygons in child nodes
    public double calculatePolygonRatio() {
        if (this.children.isEmpty()) {
            return 1; // Leaf node
        }

        int totalChildPolygons = this.children.stream().mapToInt(child -> child.polygonCount).sum();
        return totalChildPolygons == 0 ? 1 : (double) this.polygonCount / totalChildPolygons;
    }
    //HDoV search algorithm code ends

    //HDoV search algorithm code starts
    // Count the number of leaf descendants
    public int countLeafDescendants() {
        if (this.children.isEmpty()) {
            return 1; // Leaf node
        }

        return this.children.stream().mapToInt(HDoVTreeNode::countLeafDescendants).sum();
    }
    //HDoV search algorithm code ends

    //HDoV search algorithm code starts
    //private static int calculatePolygonCount(Geometry geometry) {
    public static int calculatePolygonCount(Geometry geometry) {
        int polygonCount = 0;

        if (geometry instanceof Polygon) {
            polygonCount = ((Polygon) geometry).getNumInteriorRing() + 1; // 1 for the exterior, plus any interior rings (holes)
        } else if (geometry instanceof MultiPolygon) {
            MultiPolygon multiPolygon = (MultiPolygon) geometry;
            for (int i = 0; i < multiPolygon.getNumGeometries(); i++) {
                Polygon polygon = (Polygon) multiPolygon.getGeometryN(i);
                polygonCount += polygon.getNumInteriorRing() + 1;
            }
        }
        // Add more cases if you have other types of geometries (like GeometryCollection)

        return polygonCount;
        //return 100; // Example fixed count
    }
    //HDoV search algorithm code ends

}
//HDoV Logical Structure code ends

public class VisibilityCalculator3 {
    //test code starts
    //private
    static final GeometryFactory geometryFactory = new GeometryFactory();
    //test code ends
    //private static final double DOV_THRESHOLD = 0.5; // Example threshold
    //testing with different threshold or ETA value
    private static final double DOV_THRESHOLD = 0.005; // Example threshold
    private static final String URL = "jdbc:postgresql://localhost:5432/citydb_v4";
    private static final String USER = "citydb_user";
    private static final String PASSWORD = "123123";

    public static void main(String[] args) {

        Calendar calendar = Calendar.getInstance();
        calendar.set(2000, Calendar.JANUARY, 1, 0, 0, 0);
        calendar.set(Calendar.MILLISECOND, 0);
        Date dateOfInterest = calendar.getTime();


        // Assume we have a list of JTS Geometry objects
        List<Geometry> geometryList = new ArrayList<>();
        List<Integer> buildingIds = new ArrayList<>();
        List<Date> creationDates = new ArrayList<>();
        // Vertical Storage Scheme code starts
        VPageIndex vPageIndex = new VPageIndex(100); // Assuming 100 nodes for example
        // Vertical Storage Scheme code ends

        // Populate the geometryList with actual geometry objects
        //code for indexing starts
        //note: Either use code 1 or 2 at a time, code 1 without indexing and code 2 is indexing
        //populateGeometryList(geometryList, buildingIds); //code 1 (This is original code)
        Map<Integer, Pair<Geometry, Date>> geometryMap = populateGeometryList(geometryList, buildingIds, creationDates); //code 2
        //code for indexing ends

        // Create an STRtree for efficient spatial queries
        STRtree strTree = new STRtree();

        for (Geometry geometry : geometryList) {
            strTree.insert(geometry.getEnvelopeInternal(), geometry);
        }
        strTree.build();



        // Define a viewpoint coordinate
//        Coordinate viewpoint = new Coordinate(-10000, 10000, 10000);
        double angle = Math.toRadians(45);


        // Compute the center point of all geometries
        Coordinate center = new Coordinate(0, 0, 0);
        for (Geometry geometry : geometryList) {
            Coordinate centroid = geometry.getCentroid().getCoordinate();
            center.x += centroid.x;
            center.y += centroid.y;
            center.z += centroid.z;
        }
        center.x /= geometryList.size();
        center.y /= geometryList.size();
        center.z /= geometryList.size();

// Define a distance from the center point
        double distance = 100000;

        Coordinate viewpoint = new Coordinate(
                center.x + distance * Math.cos(angle),
                center.y + distance * Math.sin(angle),
                center.z + distance
        );

        Coordinate[] viewpoints = new Coordinate[6];
        viewpoints[0] = new Coordinate(center.x + distance, center.y, center.z); // Right
        viewpoints[1] = new Coordinate(center.x - distance, center.y, center.z); // Left
        viewpoints[2] = new Coordinate(center.x, center.y + distance, center.z); // Up
        viewpoints[3] = new Coordinate(center.x, center.y - distance, center.z); // Down
        viewpoints[4] = new Coordinate(center.x, center.y, center.z + distance); // Front
        viewpoints[5] = new Coordinate(center.x, center.y, center.z - distance); // Back

// Create six bounding boxes
        Envelope[] boundingBoxes = new Envelope[6];
        for (int i = 0; i < 6; i++) {
            boundingBoxes[i] = new Envelope(new Coordinate(center.x - distance, center.y - distance, center.z - distance),
                    new Coordinate(center.x + distance, center.y + distance, center.z + distance));
        }


        String[] angleNames = {"Right", "Left", "Up", "Down", "Front", "Back"};
        for (int viewpointsCount = 0; viewpointsCount < viewpoints.length; viewpointsCount++) {
            Coordinate currentViewpoint = viewpoints[viewpointsCount];
            Envelope currentBoundingBox = boundingBoxes[viewpointsCount];
            String angleName = angleNames[viewpointsCount];
            System.out.println("=====================================================================");
            System.out.println("Viewpoint: " + angleName);

            // Calculate the visibility for each geometry object within the bounding box
            System.out.println("SRTree");
            int visibleCountSRTree = 0;
            int intersectCountSRTree = 0;
            int totalSRTreeCount = 0;
            long startTime = System.nanoTime();

            //HDoV Logical Structure code starts
            List<Pair<Integer, Double>> visibilityPairsSRTree = new ArrayList<>();
            //List<Double> visibilityValuesSRTree = new ArrayList<>();
            //HDoV Logical Structure code ends

            for (int srtreeCount = 0; srtreeCount < geometryList.size(); srtreeCount++) {
                //code for indexing starts
                //note: Either use code 1 or 2 at a time, code 1 without indexing and code 2 is indexing
                //Geometry geometry = geometryList.get(srtreeCount); //code 1 (This is original code)
                int buildingId = buildingIds.get(srtreeCount); // Get the building ID(This is original code)
                Pair<Geometry, Date> geometryDatePair = geometryMap.get(buildingId);
                Geometry geometry = geometryDatePair.getLeft();
                //code for indexing ends

                totalSRTreeCount++;
                if (currentBoundingBox.intersects(geometry.getEnvelopeInternal())) {
                    intersectCountSRTree++;
                    double visibility = calculateVisibility(geometry, currentViewpoint, strTree,2);

                    if (visibility > 0) {
                        visibleCountSRTree++;
                        //HDoV Logical Structure code starts
                        visibilityPairsSRTree.add(new ImmutablePair<>(buildingId, visibility));
                        //HDoV Logical Structure code ends
                        System.out.println("DoV of building #" + buildingId + ": " + visibility); // Print building ID and visibility
                    }
                }
            }
            //HDoV Logical Structure code starts
            // Construct HDoV-tree for SRTree after the loop
            List<HDoVTreeNode> leafNodesSRTree = new ArrayList<>();
            //List<HDoVTreeNode> leafNodesSRTree = createLeafNodes(visibilityPairsSRTree, geometryList, viewpoints);
            //HDoV search algorithm code starts
            //List<Geometry> searchResults = new ArrayList<>();
            List<GeometryWithDate> searchResultsSRTree = new ArrayList<>();
            //HDoV search algorithm code ends
            // Vertical Storage Scheme code starts
            //test code starts
            Map<Integer, Double> visibilityMap = new HashMap<>();
            //test code ends
            // Vertical Storage Scheme code ends
            for (Pair<Integer, Double> pair : visibilityPairsSRTree) {
                int buildingId = pair.getKey();
                double visibility = pair.getValue();
                //code for indexing starts
                //note: Either use code 1 or 2 at a time, code 1 without indexing and code 2 is indexing
                //Geometry geometry = null; // code 1(This is original code)
                Pair<Geometry, Date> geometryDatePair = geometryMap.get(buildingId);
                Geometry geometry = geometryDatePair.getLeft(); // code 2
                Date creationDate = geometryDatePair.getRight();
                //code for indexing ends
                // Vertical Storage Scheme code starts
                //test code starts
                visibilityMap.put(buildingId, visibility);
                //test code ends
                // Vertical Storage Scheme code ends

                // Find the geometry corresponding to the buildingId
                for (int i = 0; i < buildingIds.size(); i++) {
                    if (buildingIds.get(i) == buildingId) {
                        geometry = geometryList.get(i);
                        break;
                    }
                }

                if (geometry != null) {
                    double lod = calculateLoD(geometry, currentViewpoint, visibility);
                    //HDoV search algorithm code starts
                    //int polygonCount = calculatePolygonCount(geometry); // Calculate the polygon count for this geometry
                    int polygonCount = HDoVTreeNode.calculatePolygonCount(geometry);
                    //HDoV search algorithm code ends
                    HDoVTreeNode leafNode = new HDoVTreeNode(visibility, geometry, lod, polygonCount, creationDate);//parameter polygonCount is added for HDoV search algorithm
                    leafNodesSRTree.add(leafNode);
                }
            }
            // Create internal nodes and build the HDoV-tree
            // This part will be more complex in a real scenario where you need to define how to group leaf nodes into internal nodes
            HDoVTreeNode rootNodeSRTree = HDoVTreeNode.buildHDoVTree(leafNodesSRTree);
            //HDoVTreeNode rootNodeSRTree = new HDoVTreeNode(leafNodesSRTree);
            if (rootNodeSRTree != null) {
                rootNodeSRTree.printNode();
                System.out.println("Traversing HDoV-tree for viewpoint: " + angleNames[viewpointsCount]);
                rootNodeSRTree.traverseAndPrint(viewpoints[viewpointsCount]);
//            rootNodeSRTree.printNode();
//            System.out.println("Traversing HDoV-tree for viewpoint: " + angleName);
//            rootNodeSRTree.traverseAndPrint(currentViewpoint);
                //HDoV search algorithm code starts
                // Now set the queryPoint for the HDoV tree search
                Coordinate queryPoint = viewpoints[viewpointsCount]; // Use the current viewpoint as the query point
                //rootNode.searchHDoVTree(queryPoint, searchResults, DOV_THRESHOLD);
                long start_searchHDoVTree = System.nanoTime();
                rootNodeSRTree.searchHDoVTree(queryPoint, searchResultsSRTree, DOV_THRESHOLD);
                // Print search results
                System.out.println("Search Results for HDoV tree: ");
                for (GeometryWithDate result : searchResultsSRTree) {
                    if(result.getDate().before(dateOfInterest)) {
                        System.out.println("Geometry: " + result.getGeometry().toText() + "\nConstruction Date: " + result.getDate());
                    }
                }
                long end_searchHDoVTree = System.nanoTime();
                long duration_searchHDoVTree = end_searchHDoVTree - start_searchHDoVTree;
                double durationSeconds_searchHDoVTree = duration_searchHDoVTree / 1_000_000_000.0;
                System.out.println("Time taken to search HDoVTree: " + durationSeconds_searchHDoVTree + " seconds");
                // Vertical Storage Scheme code starts
                //test code starts
                setupVPageIndex(vPageIndex, leafNodesSRTree, visibilityMap);
                //test code ends
                // Vertical Storage Scheme code ends
                // Vertical Storage Scheme code starts
                long start_vertical_storage = System.nanoTime();
                rootNodeSRTree.searchHDoVTree(queryPoint, searchResultsSRTree, DOV_THRESHOLD, vPageIndex);
                // Print search results
                System.out.println("Search Results of vertical storage): ");
                for (GeometryWithDate result : searchResultsSRTree) {
                    System.out.println("Search Result: " + result.getGeometry().toText());
                }
                long end_vertical_storage = System.nanoTime();
                long duration_vertical_storage = end_vertical_storage - start_vertical_storage;
                double durationSeconds_vertical_storage = duration_vertical_storage / 1_000_000_000.0;
                System.out.println("Time taken to run vertical storage search: " + durationSeconds_vertical_storage + " seconds");
                long verticalStorageSize = calculateVerticalStorageSize(vPageIndex);
                System.out.println("Vertical Storage Size: " + verticalStorageSize + " bytes");
                // Vertical Storage Scheme code ends
            }
            //HDoV search algorithm code ends
            //HDoV Logical Structure code ends
            long endTime = System.nanoTime();
            long duration = endTime - startTime;
            double durationSeconds = duration / 1_000_000_000.0;
            System.out.println("Time taken to run the section: " + durationSeconds + " seconds");
            System.out.println("total: " + totalSRTreeCount);
            System.out.println("in view: " + intersectCountSRTree);
            System.out.println("visible: " + visibleCountSRTree);

            // Brute force
            System.out.println("=====================================================================");
            System.out.println("Brute Force");
            int visibleCountBruteForce = 0;
            int intersectCountBruteForce = 0;
            int totalBruteForceCount = 0;
            long startTimeBruteForce = System.nanoTime();

            //HDoV Logical Structure code starts
            List<Pair<Integer, Double>> visibilityPairsBruteForce = new ArrayList<>();
            //List<Double> visibilityValuesBruteForce = new ArrayList<>();
            //HDoV Logical Structure code ends

            for (int bruteforceCount = 0; bruteforceCount < geometryList.size(); bruteforceCount++) {
                //code for indexing starts
                //note: Either use code 1 or 2 at a time, code 1 without indexing and code 2 is indexing
                //Geometry geometry = geometryList.get(bruteforceCount); //code 1 (This is original code)
                int buildingId = buildingIds.get(bruteforceCount); // Get the building ID(This is original code)
                Pair<Geometry, Date> geometryDatePair = geometryMap.get(buildingId);
                Geometry geometry = geometryDatePair.getLeft(); //code 2
                //code for indexing ends

                totalBruteForceCount++;
                if (currentBoundingBox.intersects(geometry.getEnvelopeInternal())) {
                    intersectCountBruteForce++;
                    double visibility = calculateVisibilityBruteForce(geometry, currentViewpoint, geometryList, 2);

                    if (visibility > 0) {
                        visibleCountBruteForce++;
                        //HDoV Logical Structure code starts
                        visibilityPairsBruteForce.add(new ImmutablePair<>(buildingId, visibility));
                        //HDoV Logical Structure code ends
                        System.out.println("DoV of building #" + buildingId + ": " + visibility); // Print building ID and visibility
                    }
                }
            }
            //HDoV Logical Structure code starts
            // Construct HDoV-tree for Brute Force after the loop
            List<HDoVTreeNode> leafNodesBruteForce = new ArrayList<>();
            //List<HDoVTreeNode> leafNodesBruteForce = createLeafNodes(visibilityPairsBruteForce, geometryList, viewpoints);
            //HDoV search algorithm code starts
            List<GeometryWithDate> searchResultsBruteForce = new ArrayList<>();
            //HDoV search algorithm code ends
            // Vertical Storage Scheme code starts
            //test code starts
            visibilityMap = new HashMap<>();
            //test code ends
            // Vertical Storage Scheme code ends
            for (Pair<Integer, Double> pair : visibilityPairsBruteForce) {
                int buildingId = pair.getKey();
                double visibility = pair.getValue();
                //Geometry geometry = null;// code 1(This is original code)
                Pair<Geometry, Date> geometryDatePair = geometryMap.get(buildingId);
                Geometry geometry = geometryDatePair.getLeft(); // code 2
                Date creationDate = geometryDatePair.getRight();
                // Vertical Storage Scheme code starts
                //test code starts
                visibilityMap.put(buildingId, visibility);
                //test code ends
                // Vertical Storage Scheme code ends

                // Find the geometry corresponding to the buildingId
                for (int i = 0; i < buildingIds.size(); i++) {
                    if (buildingIds.get(i) == buildingId) {
                        geometry = geometryList.get(i);
                        break;
                    }
                }

                if (geometry != null) {
                    double lod = calculateLoD(geometry, currentViewpoint, visibility);
                    //HDoV search algorithm code starts
                    //int polygonCount = calculatePolygonCount(geometry); // Calculate the polygon count for this geometry
                    int polygonCount = HDoVTreeNode.calculatePolygonCount(geometry);
                    //HDoV search algorithm code ends
                    HDoVTreeNode leafNode = new HDoVTreeNode(visibility, geometry, lod, polygonCount, creationDate);//parameter polygonCount is added for HDoV search algorithm
                    leafNodesBruteForce.add(leafNode);
                }
            }
            // Create internal nodes and build the HDoV-tree
            // This part will be more complex in a real scenario where you need to define how to group leaf nodes into internal nodes
            HDoVTreeNode rootNodeBruteForce = HDoVTreeNode.buildHDoVTree(leafNodesBruteForce);
            //HDoVTreeNode rootNodeBruteForce = new HDoVTreeNode(leafNodesBruteForce);
            if (rootNodeBruteForce != null) {
                rootNodeBruteForce.printNode();
                System.out.println("Traversing HDoV-tree for viewpoint (Brute Force): " + angleNames[viewpointsCount]);
                rootNodeBruteForce.traverseAndPrint(viewpoints[viewpointsCount]);
                //HDoV search algorithm code starts
                // Now set the queryPoint for the HDoV tree search
                Coordinate queryPoint = viewpoints[viewpointsCount]; // Use the current viewpoint as the query point
                //rootNode.searchHDoVTree(queryPoint, searchResults, DOV_THRESHOLD);
                long start_searchHDoVTree = System.nanoTime();
                rootNodeBruteForce.searchHDoVTree(queryPoint, searchResultsBruteForce, DOV_THRESHOLD);
                // Print search results
                for (GeometryWithDate result : searchResultsBruteForce) {
                    System.out.println("Search Result: " + result.getGeometry().toText());
                }
                long end_searchHDoVTree = System.nanoTime();
                long duration_searchHDoVTree = end_searchHDoVTree - start_searchHDoVTree;
                double durationSeconds_searchHDoVTree = duration_searchHDoVTree / 1_000_000_000.0;
                System.out.println("Time taken to search HDoVTree: " + durationSeconds_searchHDoVTree + " seconds");

                // Vertical Storage Scheme code starts
                //test code starts
                setupVPageIndex(vPageIndex, leafNodesSRTree, visibilityMap);
                //test code ends
                // Vertical Storage Scheme code ends
                // Vertical Storage Scheme code starts
                long start_vertical_storage = System.nanoTime();
                rootNodeBruteForce.searchHDoVTree(queryPoint, searchResultsBruteForce, DOV_THRESHOLD, vPageIndex);
                // Print search results
                System.out.println("Search Results of vertical storage): ");
                for (GeometryWithDate result : searchResultsBruteForce) {
                    System.out.println("Search Result: " + result.getGeometry().toText());
                }
                long end_vertical_storage = System.nanoTime();
                long duration_vertical_storage = end_vertical_storage - start_vertical_storage;
                double durationSeconds_vertical_storage = duration_vertical_storage / 1_000_000_000.0;
                System.out.println("Time taken to run vertical storage search: " + durationSeconds_vertical_storage + " seconds");
                long verticalStorageSize = calculateVerticalStorageSize(vPageIndex);
                System.out.println("Vertical Storage Size: " + verticalStorageSize + " bytes");
                // Vertical Storage Scheme code ends
            }
            //HDoV search algorithm code ends
            //HDoV Logical Structure code ends
            long endTimeBruteForce = System.nanoTime();
            long durationBruteForce = endTimeBruteForce - startTimeBruteForce;
            double durationSecondsBruteForce = durationBruteForce / 1_000_000_000.0;
            System.out.println("Time taken to run the section: " + durationSecondsBruteForce + " seconds");
            System.out.println("total: " + totalBruteForceCount);
            System.out.println("in view: " + intersectCountBruteForce);
            System.out.println("visible: " + visibleCountBruteForce);

        }
    }


    //code for indexing starts

//note: Either use code 1 or 2 at a time, code 1 without indexing and code 2 is indexing

    //code 2
    private static Map<Integer, Pair<Geometry, Date>> populateGeometryList(
            List<Geometry> geometryList,
            List<Integer> buildingIds,
            List<Date> constructionDates) {

        String sql = "SELECT b.id, ST_AsText(ST_ConvexHull((ST_Dump(ms.solid_geometry)).geom)) as wkt, b.year_of_construction " +
                "FROM building b " +
                "JOIN surface_geometry ms ON b.lod2_solid_id = ms.id;";

        Map<Integer, List<Geometry>> buildingGeometries = new HashMap<>();
        Map<Integer, Date> buildingConstructionDates = new HashMap<>(); // Changed to store Date objects

        try (Connection conn = DriverManager.getConnection(URL, USER, PASSWORD);
             PreparedStatement pstmt = conn.prepareStatement(sql);
             ResultSet rs = pstmt.executeQuery()) {

            WKTReader wktReader = new WKTReader(new GeometryFactory());

            while (rs.next()) {
                int buildingId = rs.getInt("id");
                String wkt = rs.getString("wkt");
                Date constructionDate = rs.getDate("year_of_construction"); // Get the construction date

                try {
                    Geometry geometry = wktReader.read(wkt);
                    buildingGeometries.computeIfAbsent(buildingId, id -> new ArrayList<>()).add(geometry);
                    buildingConstructionDates.put(buildingId, constructionDate); // Store the construction date directly
                } catch (ParseException e) {
                    System.err.println("Error parsing WKT: " + e.getMessage());
                }
            }
        } catch (SQLException e) {
            System.err.println("Error fetching data from the database: " + e.getMessage());
        }

        Map<Integer, Pair<Geometry, Date>> geometryMap = new HashMap<>();

        for (Map.Entry<Integer, List<Geometry>> entry : buildingGeometries.entrySet()) {
            int buildingId = entry.getKey();
            List<Geometry> buildingGeometry = entry.getValue();
            Geometry combinedGeometry = UnaryUnionOp.union(buildingGeometry);

            geometryList.add(combinedGeometry);
            buildingIds.add(buildingId);
            constructionDates.add(buildingConstructionDates.get(buildingId)); // Add the date to the list

            geometryMap.put(buildingId, new ImmutablePair<>(combinedGeometry, buildingConstructionDates.get(buildingId)));
        }

        return geometryMap;
    }

    //code 1(This is original code)
//    private static void populateGeometryList(List<Geometry> geometryList, List<Integer> buildingIds) {
//        // Your database connection details
//        String sql = "SELECT b.id, ST_AsText(ST_ConvexHull((ST_Dump(ms.solid_geometry)).geom)) as wkt " +
//                "FROM building b " +
//                "JOIN surface_geometry ms ON b.lod2_solid_id = ms.id;";
//
//        Map<Integer, List<Geometry>> buildingGeometries = new HashMap<>();
//
//        try (Connection conn = DriverManager.getConnection(URL, USER, PASSWORD);
//             PreparedStatement pstmt = conn.prepareStatement(sql);
//             ResultSet rs = pstmt.executeQuery()) {
//
//            WKTReader wktReader = new WKTReader(new GeometryFactory());
//
//            while (rs.next()) {
//                int buildingId = rs.getInt("id");
//                String wkt = rs.getString("wkt");
//
//                try {
//                    Geometry geometry = wktReader.read(wkt);
//                    buildingGeometries.computeIfAbsent(buildingId, id -> new ArrayList<>()).add(geometry);
//                } catch (ParseException e) {
//                    System.err.println("Error parsing WKT: " + e.getMessage());
//                }
//            }
//        } catch (SQLException e) {
//            System.err.println("Error fetching data from the database: " + e.getMessage());
//        }
//        for (Map.Entry<Integer, List<Geometry>> entry : buildingGeometries.entrySet()) {
//            int buildingId = entry.getKey();
//            List<Geometry> buildingGeometry = entry.getValue();
//            Geometry combinedGeometry = UnaryUnionOp.union(buildingGeometry);
//            geometryList.add(combinedGeometry);
//            buildingIds.add(buildingId); // Add the building ID to the buildingIds list
//        }
//    }
//code for indexing ends

    //HDoV Logical Structure code starts
//    private static double calculateLoD(Geometry geometry, Coordinate viewpoint) {
//        double distance = geometry.getCentroid().distance(new GeometryFactory().createPoint(viewpoint));
//        // Example LoD calculation based on distance
//        if (distance < 50000) {
//            return 3; // High detail
//        } else if (distance < 100000) {
//            return 2; // Medium detail
//        } else {
//            return 1; // Low detail
//        }
//    }

    private static double calculateLoD(Geometry geometry, Coordinate viewpoint, double visibility) {
        double distance = geometry.getCentroid().distance(new GeometryFactory().createPoint(viewpoint));

        if (visibility > DOV_THRESHOLD) {
            // High detail if visibility is above the threshold
            return distance < 50000 ? 3 : 2;
        } else {
            // Low detail otherwise
            return 1;
        }
    }

    //HDoV Logical Structure code ends
    private static double calculateVisibility(Geometry geometry, Coordinate viewpoint, STRtree strTree, double gridSize) {
        List<Geometry> gridCells = createGridCells(geometry, gridSize);
        int visibleCount = 0;
        int totalCount = gridCells.size();

        // Create an STRtree for the grid cells
        STRtree cellTree = new STRtree();
        for (Geometry cell : gridCells) {
            cellTree.insert(cell.getEnvelopeInternal(), cell);
        }

        for (Geometry cell : gridCells) {
            Coordinate cellCenter = cell.getCentroid().getCoordinate();
            LineString lineOfSight = new GeometryFactory().createLineString(new Coordinate[]{viewpoint, cellCenter});

            // Query the STRtree of geometries and the STRtree of grid cells
            List<Geometry> intersectingGeometries = strTree.query(lineOfSight.getEnvelopeInternal());
            List<Geometry> intersectingCells = cellTree.query(lineOfSight.getEnvelopeInternal());

            boolean occluded = false;

            // Check for intersections with the geometries
            for (Geometry intersectingGeometry : intersectingGeometries) {
                if (intersectingGeometry.equals(geometry)) {
                    continue;
                }

                if (intersectingGeometry.intersects(lineOfSight)) {
                    occluded = true;
                    break;
                }
            }

            // If the line of sight is not occluded by other geometries,
            // check for intersections with the grid cells
            if (!occluded) {
                for (Geometry intersectingCell : intersectingCells) {
                    if (intersectingCell.equals(cell)) {
                        continue;
                    }

                    if (intersectingCell.intersects(lineOfSight)) {
                        occluded = true;
                        break;
                    }
                }
            }

            if (!occluded) {
                visibleCount++;
            }
        }

        return (double) visibleCount / totalCount;
    }


    private static double calculateVisibilityBruteForce(Geometry geometry, Coordinate viewpoint, List<Geometry> geometryList, double gridSize) {
        List<Geometry> gridCells = createGridCells(geometry, gridSize);
        int visibleCount = 0;
        int totalCount = gridCells.size();

        for (Geometry cell : gridCells) {
            Coordinate cellCenter = cell.getCentroid().getCoordinate();
            LineString lineOfSight = new GeometryFactory().createLineString(new Coordinate[]{viewpoint, cellCenter});

            boolean occluded = false;

            // Check for intersections with the geometries
            for (Geometry intersectingGeometry : geometryList) {
                if (intersectingGeometry.equals(geometry)) {
                    continue;
                }

                if (intersectingGeometry.intersects(lineOfSight)) {
                    occluded = true;
                    break;
                }
            }

            // If the line of sight is not occluded by other geometries,
            // check for intersections with the grid cells
            if (!occluded) {
                for (Geometry otherCell : gridCells) {
                    if (otherCell.equals(cell)) {
                        continue;
                    }

                    if (otherCell.intersects(lineOfSight)) {
                        occluded = true;
                        break;
                    }
                }
            }

            if (!occluded) {
                visibleCount++;
            }
        }

        return (double) visibleCount / totalCount;
    }


    private static List<Geometry> createGridCells(Geometry geometry, double gridSize) {
        Envelope envelope = geometry.getEnvelopeInternal();
        GeometryFactory geometryFactory = new GeometryFactory();
        List<Geometry> gridCells = new ArrayList<>();

        for (double x = envelope.getMinX(); x < envelope.getMaxX(); x += gridSize) {
            for (double y = envelope.getMinY(); y < envelope.getMaxY(); y += gridSize) {
                Polygon cell = geometryFactory.createPolygon(new Coordinate[]{
                        new Coordinate(x, y),
                        new Coordinate(x + gridSize, y),
                        new Coordinate(x + gridSize, y + gridSize),
                        new Coordinate(x, y + gridSize),
                        new Coordinate(x, y)
                });

                if (cell.intersects(geometry)) {
                    Geometry intersection = cell.intersection(geometry);
                    gridCells.add(intersection);
                }
            }
        }

        return gridCells;
    }

//    private static List<HDoVTreeNode> createLeafNodes(List<Pair<Integer, Double>> visibilityPairs, List<Geometry> geometryList, List<Integer> buildingIds, Coordinate[] viewpoints) {
//        List<HDoVTreeNode> leafNodes = new ArrayList<>();
//        for (Pair<Integer, Double> pair : visibilityPairs) {
//            int buildingId = pair.getKey();
//            double visibility = pair.getValue();
//            Geometry geometry = null;
//
//            // Find the geometry corresponding to the buildingId
//            for (int i = 0; i < geometryList.size(); i++) {
//                if (buildingIds.get(i) == buildingId) {
//                    geometry = geometryList.get(i);
//                    break;
//                }
//            }
//
//            if (geometry != null) {
//                double lod = calculateLoD(geometry, viewpoints[0]); // Assuming first viewpoint for simplicity
//                leafNodes.add(new HDoVTreeNode(visibility, geometry, lod));
//            }
//        }
//        return leafNodes;
//    }

    // Vertical Storage Scheme code starts
    //test code starts
    private static void setupVPageIndex(VPageIndex vPageIndex, List<HDoVTreeNode> leafNodes, Map<Integer, Double> visibilityMap) {
        // Loop through each node and set up V-page-index
        int offset = 0;
        for (HDoVTreeNode node : leafNodes) {
            if (offset >= vPageIndex.getSize()) {

                //removal of code to address problem of "Reached the limit of VPageIndex. Some nodes will not be processed" starts
                //System.out.println("Warning: Reached the limit of VPageIndex. Some nodes will not be processed.");
                //break; // Break out of the loop when the limit is reached
                //removal of code to address problem of "Reached the limit of VPageIndex. Some nodes will not be processed" ends

                //Addition of code to address problem of "Reached the limit of VPageIndex. Some nodes will not be processed" starts
                vPageIndex.resize(vPageIndex.getSize() + 100); // Increase by 100, adjust as needed
                //Addition of code to address problem of "Reached the limit of VPageIndex. Some nodes will not be processed" ends
            }
            node.setVPageIndexOffset(offset);
            VPage vPage = new VPage();
            // Populate vPage with visibility data for this node
            //double visibilityValue = calculateNodeVisibility(node); // Implement this method based on your criteria
            Double visibilityValue = visibilityMap.getOrDefault(node.getId(), 0.0);
            vPage.addVisibilityData(node.getId(), visibilityValue);

            vPageIndex.setVPage(offset, vPage);
            offset++;
        }
    }
    //test code ends
    // Vertical Storage Scheme code ends

    // Vertical Storage Scheme code starts
    // Method to calculate the storage size of the vertical scheme
    private static long calculateVerticalStorageSize(VPageIndex vPageIndex) {
        long size = 0;

        // Calculate the storage size of the VPageIndex
        size += vPageIndex.getSize() * sizeOfPointer();

        // Calculate the storage size of each VPage
        for (VPage vPage : vPageIndex.getVPages()) {
            if (vPage != null) {
                size += 4; // Size of the HashMap object itself (approximation)
                size += vPage.getVisibilityData().size() * (4 + 8); // Size of each entry (integer key + double value)
            }
        }

        return size;
    }

    // Helper method to get the size of a pointer
    private static int sizeOfPointer() {
        return Integer.SIZE / Byte.SIZE; // Replace with 8 if on a 64-bit architecture
    }
    // Vertical Storage Scheme code ends

}
