package hdov;
import java.sql.*;
import java.time.LocalDate;
import java.util.concurrent.ThreadLocalRandom;

public class PopulateCreationDate {

    public static void main(String[] args) {
        String url = "jdbc:postgresql://localhost:5432/citydb_v4";
        String user = "citydb_user";
        String password = "123123";

        try (Connection conn = DriverManager.getConnection(url, user, password)) {
            System.out.println("Connected to the PostgreSQL server successfully.");

            String selectQuery = "SELECT id FROM building";
            Statement stmt = conn.createStatement();
            ResultSet rs = stmt.executeQuery(selectQuery);

            String updateQuery = "UPDATE building SET year_of_construction = ? WHERE id = ?";
            PreparedStatement updateStmt = conn.prepareStatement(updateQuery);

            while (rs.next()) {
                int id = rs.getInt("id");

                long minDay = LocalDate.of(1980, 1, 1).toEpochDay();
                long maxDay = LocalDate.of(2010, 12, 31).toEpochDay();
                long randomDay = ThreadLocalRandom.current().nextLong(minDay, maxDay + 1);
                LocalDate randomCreationDate = LocalDate.ofEpochDay(randomDay);

                updateStmt.setDate(1, java.sql.Date.valueOf(randomCreationDate));
                updateStmt.setInt(2, id);

                updateStmt.executeUpdate();
            }
            updateStmt.close();
            stmt.close();
            System.out.println("Creation dates updated successfully with unique random values.");
        } catch (SQLException e) {
            System.out.println(e.getMessage());
        }
    }
}
