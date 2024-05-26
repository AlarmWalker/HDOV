package hdov;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class VisibilityCalculator2 {
    private static final String URL = "jdbc:postgresql://localhost:5432/citydb_v4";
    private static final String USER = "citydb_user";
    private static final String PASSWORD = "123123";

    public static void main(String[] args) {
        // Assume we have a list of JTS Geometry objects 
        List<Geometry> geometryList = new ArrayList<>();
        List<Integer> buildingIds = new ArrayList<>();

        // Populate the geometryList with actual geometry objects
        populateGeometryList(geometryList, buildingIds);

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

            for (int srtreeCount = 0; srtreeCount < geometryList.size(); srtreeCount++) {
                Geometry geometry = geometryList.get(srtreeCount);
                int buildingId = buildingIds.get(srtreeCount); // Get the building ID

                totalSRTreeCount++;
                if (currentBoundingBox.intersects(geometry.getEnvelopeInternal())) {
                    intersectCountSRTree++;
                    double visibility = calculateVisibility(geometry, currentViewpoint, strTree,2);

                    if (visibility > 0) {
                        visibleCountSRTree++;
                        System.out.println("DoV of building #" + buildingId + ": " + visibility); // Print building ID and visibility
                    }
                }
            }
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

            for (int bruteforceCount = 0; bruteforceCount < geometryList.size(); bruteforceCount++) {
                Geometry geometry = geometryList.get(bruteforceCount);
                int buildingId = buildingIds.get(bruteforceCount); // Get the building ID

                totalBruteForceCount++;
                if (currentBoundingBox.intersects(geometry.getEnvelopeInternal())) {
                    intersectCountBruteForce++;
                    double visibility = calculateVisibilityBruteForce(geometry, currentViewpoint, geometryList, 2);

                    if (visibility > 0) {
                        visibleCountBruteForce++;
                        System.out.println("DoV of building #" + buildingId + ": " + visibility); // Print building ID and visibility
                    }
                }
            }
            long endTimeBruteForce = System.nanoTime();
            long durationBruteForce = endTimeBruteForce - startTimeBruteForce;
            double durationSecondsBruteForce = durationBruteForce / 1_000_000_000.0;
            System.out.println("Time taken to run the section: " + durationSecondsBruteForce + " seconds");
            System.out.println("total: " + totalBruteForceCount);
            System.out.println("in view: " + intersectCountBruteForce);
            System.out.println("visible: " + visibleCountBruteForce);

        }
    }


    private static void populateGeometryList(List<Geometry> geometryList, List<Integer> buildingIds) {
        // Your database connection details


        String sql = "SELECT b.id, ST_AsText(ST_ConvexHull((ST_Dump(ms.solid_geometry)).geom)) as wkt " +
                "FROM building b " +
                "JOIN surface_geometry ms ON b.lod2_solid_id = ms.id;";

        Map<Integer, List<Geometry>> buildingGeometries = new HashMap<>();

        try (Connection conn = DriverManager.getConnection(URL, USER, PASSWORD);
             PreparedStatement pstmt = conn.prepareStatement(sql);
             ResultSet rs = pstmt.executeQuery()) {

            WKTReader wktReader = new WKTReader(new GeometryFactory());

            while (rs.next()) {
                int buildingId = rs.getInt("id");
                String wkt = rs.getString("wkt");

                try {
                    Geometry geometry = wktReader.read(wkt);
                    buildingGeometries.computeIfAbsent(buildingId, id -> new ArrayList<>()).add(geometry);
                } catch (ParseException e) {
                    System.err.println("Error parsing WKT: " + e.getMessage());
                }
            }
        } catch (SQLException e) {
            System.err.println("Error fetching data from the database: " + e.getMessage());
        }
        for (Map.Entry<Integer, List<Geometry>> entry : buildingGeometries.entrySet()) {
            int buildingId = entry.getKey();
            List<Geometry> buildingGeometry = entry.getValue();
            Geometry combinedGeometry = UnaryUnionOp.union(buildingGeometry);
            geometryList.add(combinedGeometry);
            buildingIds.add(buildingId); // Add the building ID to the buildingIds list
        }
    }

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

}


class TemporalData {
    Geometry geometry;
    Timestamp timestamp;

    public TemporalData(Geometry geometry, Timestamp timestamp) {
        this.geometry = geometry;
        this.timestamp = timestamp;
    }
}

class BPlusTree {
    private Node root;
    private static final int MAX_ENTRIES = 5;
    private static final int MAX_CHILDREN = 4;

    public BPlusTree() {
        this.root = new LeafNode();
    }

    public void insert(double key, TemporalData data) {
        Split result = root.insert(key, data);
        if (result != null) {
            InternalNode newRoot = new InternalNode();
            newRoot.keys.add(result.key);
            newRoot.children.add(result.left);
            newRoot.children.add(result.right);
            this.root = newRoot;
        }
    }

    public List<TemporalData> rangeQuery(double startKey, double endKey) {
        return root.rangeQuery(startKey, endKey);
    }

    private abstract class Node {
        List<Double> keys = new ArrayList<>();

        abstract Split insert(double key, TemporalData data);
        abstract List<TemporalData> rangeQuery(double startKey, double endKey);
    }

    private class Split {
        double key;
        Node left;
        Node right;

        Split(double k, Node l, Node r) {
            this.key = k;
            this.left = l;
            this.right = r;
        }
    }

    private class InternalNode extends Node {
        List<Node> children = new ArrayList<>();

        @Override
        Split insert(double key, TemporalData data) {
            int loc = findLocation(key);
            Split split = children.get(loc).insert(key, data);
            if (split == null) return null;

            int idx = findLocation(split.key);
            keys.add(idx, split.key);
            children.add(idx + 1, split.right);

            if (keys.size() > MAX_CHILDREN) {
                int mid = keys.size() / 2;
                InternalNode sibling = new InternalNode();
                sibling.keys.addAll(keys.subList(mid + 1, keys.size()));
                sibling.children.addAll(children.subList(mid + 1, children.size() + 1));

                keys.subList(mid, keys.size()).clear();
                children.subList(mid + 1, children.size()).clear();

                return new Split(keys.remove(mid), this, sibling);
            }
            return null;
        }

        @Override
        List<TemporalData> rangeQuery(double startKey, double endKey) {
            List<TemporalData> result = new ArrayList<>();
            for (int i = 0; i < keys.size(); i++) {
                if (startKey < keys.get(i)) {
                    result.addAll(children.get(i).rangeQuery(startKey, endKey));
                }
                if (endKey < keys.get(i)) break;
            }
            result.addAll(children.get(children.size() - 1).rangeQuery(startKey, endKey));
            return result;
        }

        private int findLocation(double key) {
            int loc = 0;
            while (loc < keys.size() && keys.get(loc) < key) {
                loc++;
            }
            return loc;
        }
    }

    private class LeafNode extends Node {
        List<TemporalData> values = new ArrayList<>();
        LeafNode next;

        @Override
        Split insert(double key, TemporalData data) {
            int loc = findLocation(key);
            keys.add(loc, key);
            values.add(loc, data);

            if (keys.size() > MAX_ENTRIES) {
                LeafNode sibling = new LeafNode();
                int mid = (keys.size() + 1) / 2;

                sibling.keys.addAll(keys.subList(mid, keys.size()));
                sibling.values.addAll(values.subList(mid, values.size()));

                keys.subList(mid, keys.size()).clear();
                values.subList(mid, values.size()).clear();
                sibling.next = next;
                next = sibling;

                return new Split(sibling.keys.get(0), this, sibling);
            }
            return null;
        }

        @Override
        List<TemporalData> rangeQuery(double startKey, double endKey) {
            List<TemporalData> result = new ArrayList<>();
            LeafNode current = this;
            while (current != null) {
                for (int i = 0; i < current.keys.size(); i++) {
                    if (current.keys.get(i) > endKey) return result;
                    if (current.keys.get(i) >= startKey && current.keys.get(i) <= endKey) {
                        result.add(current.values.get(i));
                    }
                }
                current = current.next;
            }
            return result;
        }

        private int findLocation(double key) {
            int loc = 0;
            while (loc < keys.size() && keys.get(loc) < key) {
                loc++;
            }
            return loc;
        }
    }
}
