document.addEventListener("DOMContentLoaded", function () {
  // Get all navigation toggles
  const navToggles = document.querySelectorAll(".nav-toggle");

  navToggles.forEach(function (toggle) {
    toggle.addEventListener("click", function () {
      // Get the associated nav-links
      const navLinks = this.nextElementSibling;
      const arrow = this.querySelector(".arrow");

      // Toggle active class
      this.classList.toggle("active");
      navLinks.classList.toggle("active");

      // Rotate arrow
      if (this.classList.contains("active")) {
        arrow.style.transform = "rotate(90deg)";
      } else {
        arrow.style.transform = "rotate(0deg)";
      }
    });
  });
});
